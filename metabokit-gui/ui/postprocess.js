const $ = (id) => document.getElementById(id);

function number(id) {
  const value = Number($(id).value);
  return Number.isFinite(value) ? Math.max(0, value) : 0;
}

function format(value, digits = 2, suffix = "") {
  return Number.isFinite(value) ? `${Number(value).toFixed(digits)}${suffix}` : "—";
}

export function createPostprocessor({ invoke, saveDialog, asPath, log }) {
  const els = {
    empty: $("post-empty"),
    workspace: $("post-workspace"),
    report: $("post-report"),
    status: $("post-status"),
    preview: $("post-preview"),
    export: $("post-export"),
    reset: $("post-reset"),
    stats: $("post-stats"),
    tableWrap: $("post-table-wrap"),
    rows: $("post-rows"),
  };
  const filterIds = [
    "post-detected", "post-shape", "post-score", "post-sn", "post-peaks",
    "post-cv", "post-identified", "post-msms", "post-isf",
  ];

  let outputDir = null;
  let active = false;
  let running = false;
  let loaded = false;
  let busy = false;
  let hasPreview = false;
  let requestToken = 0;

  function options() {
    const cvText = $("post-cv").value.trim();
    const cv = cvText === "" ? null : number("post-cv");
    return {
      minimumDetectedPercent: Math.min(100, number("post-detected")),
      minimumPeakShape: Math.min(1, number("post-shape")),
      minimumScore: Math.min(1, number("post-score")),
      minimumSn: number("post-sn"),
      minimumMatchingPeaks: number("post-peaks"),
      maximumCvPercent: cv && cv > 0 ? cv : null,
      identifiedOnly: $("post-identified").checked,
      msmsOnly: $("post-msms").checked,
      removeIsf: $("post-isf").checked,
    };
  }

  function setEmpty(title, detail) {
    els.empty.querySelector("p").textContent = title;
    els.empty.querySelector("span").textContent = detail;
    els.empty.classList.remove("is-hidden");
    els.workspace.classList.add("is-hidden");
  }

  function invalidate() {
    requestToken += 1;
    hasPreview = false;
    els.export.disabled = true;
    els.status.textContent = "Filters changed · preview to update the table.";
  }

  function setBusy(on, message = "") {
    busy = on;
    if (message) els.status.textContent = message;
    [els.report, els.preview, els.reset, ...filterIds.map($)].forEach((element) => {
      element.disabled = on || running;
    });
    if (on || running) els.export.disabled = true;
  }

  async function open() {
    if (!active || !invoke || !outputDir || running) return;
    const token = ++requestToken;
    setBusy(true, "Finding completed reports…");
    try {
      const session = await invoke("postprocess_open", { outputDir });
      if (token !== requestToken) return;
      els.report.textContent = "";
      for (const report of session.reports) {
        const option = document.createElement("option");
        option.value = report;
        option.textContent = report;
        els.report.append(option);
      }
      loaded = session.reports.length > 0;
      if (!loaded) {
        setEmpty("No completed reports found", "Run the analysis first, then return here to clean its tables.");
        return;
      }
      els.empty.classList.add("is-hidden");
      els.workspace.classList.remove("is-hidden");
      setBusy(false);
      await preview();
    } catch (error) {
      if (token !== requestToken) return;
      loaded = false;
      setEmpty("Could not open the reports", String(error));
      log(String(error), "error");
    } finally {
      if (token === requestToken) setBusy(false);
    }
  }

  function render(result) {
    hasPreview = true;
    els.stats.textContent = "";
    const removed = result.totalRows - result.keptRows;
    const items = [
      ["Source rows", result.totalRows],
      ["Rows kept", result.keptRows],
      ["Rows removed", removed],
    ];
    for (const [label, value] of items) {
      const wrap = document.createElement("div");
      const dt = document.createElement("dt");
      const dd = document.createElement("dd");
      dt.textContent = label;
      dd.textContent = String(value);
      wrap.append(dt, dd);
      els.stats.append(wrap);
    }
    if (result.removedBy.length) {
      const detail = document.createElement("p");
      detail.className = "mono micro muted post-removals";
      detail.textContent = result.removedBy
        .map((item) => `${item.criterion}: ${item.rows}`)
        .join(" · ");
      els.stats.append(detail);
    }

    els.rows.textContent = "";
    for (const row of result.rows) {
      const tr = document.createElement("tr");
      const values = [
        row.group,
        row.name || "Unidentified",
        row.adduct,
        format(row.mz, 4),
        format(row.detectedPercent, 0, "%"),
        format(row.peakShape, 2),
        format(row.score, 2),
        format(row.sn, 1),
        format(row.cvPercent, 1, "%"),
      ];
      values.forEach((value, index) => {
        const td = document.createElement("td");
        td.textContent = value;
        if (index === 1) td.title = value;
        tr.append(td);
      });
      els.rows.append(tr);
    }
    els.tableWrap.classList.toggle("is-hidden", result.rows.length === 0);
    els.status.textContent = result.previewTruncated
      ? `Showing the first ${result.rows.length} of ${result.keptRows} retained rows.`
      : `${result.keptRows} retained row${result.keptRows === 1 ? "" : "s"}.`;
  }

  async function preview() {
    if (!loaded || busy || running) return;
    const token = ++requestToken;
    setBusy(true, "Applying filters…");
    try {
      const result = await invoke("postprocess_preview", {
        outputDir,
        report: els.report.value,
        options: options(),
      });
      if (token !== requestToken) return;
      render(result);
      els.export.disabled = false;
    } catch (error) {
      if (token !== requestToken) return;
      els.status.textContent = "Could not preview this report.";
      log(String(error), "error");
    } finally {
      if (token === requestToken) {
        setBusy(false);
        els.export.disabled = running || !hasPreview;
      }
    }
  }

  async function exportCleaned() {
    if (!loaded || busy || running) return;
    const stem = els.report.value.replace(/\.csv$/i, "");
    const destination = asPath(await saveDialog({
      title: "Export cleaned report",
      defaultPath: `${stem}_clean.csv`,
      filters: [{ name: "CSV table", extensions: ["csv"] }],
    }));
    if (!destination) return;
    setBusy(true, "Writing cleaned report…");
    try {
      const result = await invoke("postprocess_export", {
        outputDir,
        report: els.report.value,
        destination,
        options: options(),
      });
      render(result);
      els.status.textContent = `Exported ${result.keptRows} rows to ${destination}.`;
      log(`exported cleaned report to ${destination}`, "done");
    } catch (error) {
      els.status.textContent = "Could not export the cleaned report.";
      log(String(error), "error");
    } finally {
      setBusy(false);
      els.export.disabled = false;
    }
  }

  function reset() {
    $("post-detected").value = "0";
    $("post-shape").value = "0";
    $("post-score").value = "0";
    $("post-sn").value = "0";
    $("post-peaks").value = "0";
    $("post-cv").value = "";
    $("post-identified").checked = false;
    $("post-msms").checked = false;
    $("post-isf").checked = true;
    preview();
  }

  filterIds.forEach((id) => $(id).addEventListener("change", invalidate));
  els.report.addEventListener("change", preview);
  els.preview.addEventListener("click", preview);
  els.export.addEventListener("click", exportCleaned);
  els.reset.addEventListener("click", reset);

  return {
    async setOutput(path, { reload = false } = {}) {
      if (!reload && path === outputDir) return;
      outputDir = path || null;
      requestToken += 1;
      loaded = false;
      hasPreview = false;
      els.rows.textContent = "";
      els.stats.textContent = "";
      els.export.disabled = true;
      if (!outputDir) {
        setEmpty("No completed reports found", "Choose a processed dataset first.");
      } else if (active) {
        await open();
      }
    },
    async activate() {
      active = true;
      if (outputDir && !loaded) await open();
    },
    setRunning(on) {
      running = on;
      setBusy(busy, on ? "Analysis is writing the report…" : "Ready to preview.");
      if (!on && active && outputDir) open();
    },
  };
}
