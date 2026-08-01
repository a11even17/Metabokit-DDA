import { createVisualizer } from "./visualizer.js";

/* MetaboKit DDA — interface logic.
 *
 * No framework and no build step. The flow is: pick a folder, the Rust side
 * scans it and returns both a ready-to-run configuration and an auditable list
 * of what it concluded; the UI shows those conclusions and only asks for input
 * where the scan came up short. Everything the old parameter tabs exposed still
 * exists, under Advanced.
 */

const TAURI = window.__TAURI__;
const invoke = TAURI ? TAURI.core.invoke : null;
const listen = TAURI ? TAURI.event.listen : null;

const LAST_DATASET_KEY = "mk-last-dataset";

const STAGES = [
  ["preparing", "Preparing"],
  ["library", "Reading libraries"],
  ["processing", "Processing samples"],
  ["aligning", "Aligning samples"],
  ["reporting", "Writing reports"],
  ["gap-filling", "Gap filling"],
];

let params = null;
let scan = null;
let running = false;
let lastOutcome = null;
let visualizer = null;

const $ = (id) => document.getElementById(id);
const $$ = (sel) => Array.from(document.querySelectorAll(sel));

/* ------------------------------------------------------------------ utils */

function bytes(n) {
  if (!n) return "—";
  const units = ["B", "KB", "MB", "GB", "TB"];
  let i = 0;
  let v = n;
  while (v >= 1024 && i < units.length - 1) {
    v /= 1024;
    i += 1;
  }
  return `${v < 10 && i > 0 ? v.toFixed(1) : Math.round(v)} ${units[i]}`;
}

function baseName(p) {
  const parts = String(p).split(/[\\/]/);
  return parts[parts.length - 1] || p;
}

/** Normalise the several shapes a Tauri dialog can return into a path string. */
function asPath(x) {
  if (!x) return null;
  if (typeof x === "string") return x;
  return x.path || null;
}

function clear(el) {
  el.textContent = "";
}

/** Append a row of spans. Text only — nothing here ever parses HTML. */
function row(list, cells, className) {
  const li = document.createElement("li");
  if (className) li.className = className;
  cells.forEach(([cls, text, title]) => {
    const span = document.createElement("span");
    if (cls) span.className = cls;
    span.textContent = text;
    if (title) span.title = title;
    li.append(span);
  });
  list.append(li);
  return li;
}

function log(message, level = "info") {
  const list = $("log");
  const li = document.createElement("li");
  li.className = level;
  const lv = document.createElement("span");
  lv.className = "lv";
  lv.textContent = level === "info" ? "" : level;
  const body = document.createElement("span");
  body.textContent = message;
  li.append(lv, body);
  list.append(li);
  const atBottom = list.scrollHeight - list.scrollTop - list.clientHeight < 40;
  if (atBottom) list.scrollTop = list.scrollHeight;
  while (list.childElementCount > 2000) list.removeChild(list.firstChild);
}

/* --------------------------------------------------------------- dialogs */

async function dialogOpen(options) {
  if (!invoke) return null;
  try {
    return await invoke("plugin:dialog|open", { options });
  } catch (err) {
    log(`file dialog unavailable: ${err}`, "error");
    return null;
  }
}

async function dialogSave(options) {
  if (!invoke) return null;
  try {
    return await invoke("plugin:dialog|save", { options });
  } catch (err) {
    log(`file dialog unavailable: ${err}`, "error");
    return null;
  }
}

async function openOutputDirectory(path) {
  if (!invoke || !path) return;
  try {
    await invoke("open_output_directory", { path });
  } catch (err) {
    log(String(err), "error");
  }
}

function storedDataset() {
  try {
    return localStorage.getItem(LAST_DATASET_KEY);
  } catch (err) {
    return null;
  }
}

function rememberDataset(path) {
  try {
    if (path) localStorage.setItem(LAST_DATASET_KEY, path);
    else localStorage.removeItem(LAST_DATASET_KEY);
  } catch (err) {
    /* Storage may be disabled; selection still works for this session. */
  }
}

/* ------------------------------------------------------------- discovery */

async function chooseDataset() {
  const dir = asPath(
    await dialogOpen({ directory: true, title: "Choose dataset folder" })
  );
  if (dir) await runScan(dir);
}

async function runScan(dir, { restored = false } = {}) {
  $("rail-stage").textContent = "Scanning…";
  try {
    scan = await invoke("scan_dataset", { path: dir });
  } catch (err) {
    $("rail-stage").textContent = "Scan failed";
    log(String(err), "error");
    return false;
  }
  rememberDataset(scan.root);
  params = scan.params;
  resetRunView();
  renderDataset();
  render();
  await refreshSystem();
  await revalidate();
  if (scan.previousRun) {
    renderPreviousRun(scan.previousRun);
    $("rail-stage").textContent = "Already processed";
  } else {
    $("rail-stage").textContent = scan.ready
      ? `${scan.samples.length} samples ready`
      : "Needs attention";
  }
  if (visualizer) {
    await visualizer.setOutput(params.outputDir, { reload: true });
  }
  log(restored ? `restored ${scan.root}` : `scanned ${scan.root}`);
  return true;
}

function renderDataset() {
  const has = Boolean(scan);
  $("pickzone").classList.toggle("is-hidden", has);
  $("dataset-detail").classList.toggle("is-hidden", !has);
  if (!has) return;

  $("ds-root").textContent = scan.root;
  $("ds-output").textContent = scan.outputDir;
  $("ds-total").textContent = `${scan.samples.length} file${
    scan.samples.length === 1 ? "" : "s"
  } · ${bytes(scan.totalBytes)}`;

  const findings = $("findings");
  clear(findings);
  const wording = { ok: "found", warn: "check", blocked: "missing" };
  scan.notes.forEach((n) => {
    row(
      findings,
      [
        ["lvl", wording[n.level] || n.level],
        ["topic", n.topic],
        ["msg", n.message],
      ],
      n.level
    );
  });

  const samples = $("samples");
  clear(samples);
  $("samples-title").textContent =
    scan.samples.length === 1 ? "1 sample" : `${scan.samples.length} samples`;
  scan.samples.forEach((s) => {
    row(samples, [
      ["fname", s.name, s.path],
      [
        "fsize",
        s.subfolder ? `${s.subfolder}  ·  ${bytes(s.bytes)}` : bytes(s.bytes),
      ],
    ]);
  });

  renderLibraries();
  refreshFixups();
}

function renderLibraries() {
  const list = $("libraries");
  clear(list);
  const libs = (params && params.libraries) || [];
  $("libraries-empty").classList.toggle("is-hidden", libs.length > 0);

  libs.forEach((lib) => {
    const li = row(list, [
      [
        "fname",
        lib.kind === "builtin" ? lib.value : baseName(lib.value),
        lib.value,
      ],
      ["fsize", lib.kind],
    ]);
    const drop = document.createElement("button");
    drop.className = "drop";
    drop.textContent = "×";
    drop.title = "Remove";
    drop.addEventListener("click", () => {
      params.libraries = params.libraries.filter(
        (l) => !(l.kind === lib.kind && l.value === lib.value)
      );
      renderLibraries();
      refreshFixups();
      revalidate();
    });
    li.append(drop);
  });
}

async function addLibraries(kind) {
  const picked = await dialogOpen({
    multiple: true,
    title: kind === "csv" ? "Add CSV library" : "Add MSP library",
    filters:
      kind === "csv"
        ? [{ name: "CSV library", extensions: ["csv"] }]
        : [{ name: "MSP library", extensions: ["msp", "txt"] }],
  });
  if (!picked) return;
  const paths = (Array.isArray(picked) ? picked : [picked])
    .map(asPath)
    .filter(Boolean);
  if (!paths.length) return;
  params = await invoke("add_libraries", { params: collect(), paths, kind });
  renderLibraries();
  refreshFixups();
  await revalidate();
}

async function chooseLibsDir() {
  const dir = asPath(
    await dialogOpen({ directory: true, title: "Locate libs/ folder" })
  );
  if (!dir) return;
  const next = collect();
  next.libsDir = dir;
  const [updated, detected] = await invoke("relink_libraries", { params: next });
  params = updated;
  $("libsDir").value = dir;
  renderLibraries();
  refreshFixups();
  await revalidate();
  log(
    detected.length
      ? `${detected.length} librar${
          detected.length === 1 ? "y" : "ies"
        } linked from ${dir}`
      : `no usable libraries in ${dir}`,
    detected.length ? "info" : "warn"
  );
}

/** Hide the fix-up strip once the thing it was offering to fix is resolved. */
function refreshFixups() {
  const stillBlocked =
    !params || !params.libraries || params.libraries.length === 0;
  $("fixups").classList.toggle("is-hidden", !stillBlocked);
}

/* ------------------------------------------------------------ param bind */

function collect() {
  if (!params) return params;
  $$("[data-param]").forEach((el) => {
    const key = el.dataset.param;
    if (el.classList.contains("seg")) return;
    if (el.type === "checkbox") {
      params[key] = el.checked;
    } else if (el.type === "number") {
      const v = parseFloat(el.value);
      params[key] = Number.isFinite(v) ? v : 0;
    } else {
      params[key] = el.value.trim();
    }
  });

  params.peakWidth = [
    parseFloat($("peakWidthMin").value) || 0,
    parseFloat($("peakWidthMax").value) || 0,
  ];
  params.excludeAdducts = $("excludeAdducts")
    .value.split(",")
    .map((x) => x.trim())
    .filter(Boolean);
  params.intensityCutoff = $("intensityCutoffOn").checked
    ? parseFloat($("intensityCutoff").value) || 0
    : null;
  params.libsDir = $("libsDir").value.trim() || null;
  params.minPeaks = Math.max(0, Math.min(255, Math.round(params.minPeaks || 0)));
  params.threads = Math.max(0, Math.round(params.threads || 0));
  params.maxFilesInFlight = Math.max(0, Math.round(params.maxFilesInFlight || 0));
  params.matchNFragments = Math.max(1, Math.round(params.matchNFragments || 1));
  return params;
}

function render() {
  if (!params) return;
  $$("[data-param]").forEach((el) => {
    const key = el.dataset.param;
    const value = params[key];
    if (el.classList.contains("seg")) {
      Array.from(el.children).forEach((b) =>
        b.classList.toggle("is-on", b.dataset.value === value)
      );
    } else if (el.type === "checkbox") {
      el.checked = Boolean(value);
    } else if (value !== null && value !== undefined) {
      el.value = value;
    } else {
      el.value = "";
    }
  });

  $("peakWidthMin").value = params.peakWidth[0];
  $("peakWidthMax").value = params.peakWidth[1];
  $("excludeAdducts").value = (params.excludeAdducts || []).join(", ");
  const hasCutoff =
    params.intensityCutoff !== null && params.intensityCutoff !== undefined;
  $("intensityCutoffOn").checked = hasCutoff;
  $("intensityCutoff").disabled = !hasCutoff;
  $("intensityCutoff").value = hasCutoff ? params.intensityCutoff : "";
  $("libsDir").value = params.libsDir || "";
  renderLibraries();
}

async function refreshSystem() {
  if (!invoke || !params) return;
  try {
    const info = await invoke("system_info", { params: collect() });
    $("about-version").textContent = info.version;
    $("about-os").textContent = info.os;
    $("about-threads").textContent = String(info.availableThreads);
    $("masthead-tag").textContent = `DDA · v${info.version}`;
    $("threads-hint").textContent = `0 uses all ${info.availableThreads} cores.`;
  } catch (err) {
    /* leave the placeholders */
  }
}

/* ------------------------------------------------------------ validation */

async function revalidate() {
  if (!invoke || !params) {
    $("btn-run").disabled = true;
    return [];
  }
  let problems = [];
  try {
    problems = await invoke("validate_params", { params: collect() });
  } catch (err) {
    problems = [{ field: "params", message: String(err), fatal: true }];
  }
  const box = $("issues");
  clear(box);
  box.classList.toggle("is-hidden", problems.length === 0);
  problems.forEach((p) => {
    const el = document.createElement("p");
    const tag = document.createElement("span");
    tag.className = `tag ${p.fatal ? "bad" : "warn"}`;
    tag.textContent = p.fatal ? "blocking" : "check";
    const msg = document.createElement("span");
    msg.textContent = p.message;
    el.append(tag, msg);
    box.append(el);
  });

  const blocked = problems.some((p) => p.fatal);
  $("btn-run").disabled = blocked || running;
  return problems;
}

/* ----------------------------------------------------------------- views */

function showView(name) {
  $$(".view").forEach((v) => v.classList.remove("is-active"));
  const view = $(`view-${name}`);
  if (view) view.classList.add("is-active");
  $$(".nav-item").forEach((b) =>
    b.classList.toggle("is-active", b.dataset.view === name)
  );
  document.querySelector(".stage").scrollTop = 0;
  if (name === "visualizer" && visualizer) {
    visualizer.activate().catch((err) => log(String(err), "error"));
  }
}

function renderStages(current) {
  const list = $("stages");
  clear(list);
  const at = current === "complete"
    ? STAGES.length
    : STAGES.findIndex(([id]) => id === current);
  STAGES.forEach(([id, label], i) => {
    const li = document.createElement("li");
    if (at >= 0 && i < at) li.classList.add("is-done");
    if (i === at) li.classList.add("is-now");
    const ord = document.createElement("span");
    ord.className = "st-ord";
    ord.textContent = String(i + 1).padStart(2, "0");
    const name = document.createElement("span");
    name.className = "st-name";
    name.textContent = label;
    li.append(ord, name);
    list.append(li);
  });
  const done = current === "complete" ? 1 : at < 0 ? 0 : (at + 1) / STAGES.length;
  $("rail-meter").style.width = `${Math.round(done * 100)}%`;
}

function renderStats(entries) {
  const stats = $("stats");
  clear(stats);
  entries.forEach(([k, v]) => {
    const wrap = document.createElement("div");
    const dt = document.createElement("dt");
    dt.textContent = k;
    const dd = document.createElement("dd");
    dd.textContent = String(v);
    wrap.append(dt, dd);
    stats.append(wrap);
  });
}

function resetRunView() {
  lastOutcome = null;
  $("result").classList.add("is-hidden");
  $("run-title").textContent = "Execution";
  $("run-lede").textContent = "Nothing has run yet for this dataset.";
  $("btn-run").textContent = "Run analysis";
  $("viz-status").textContent =
    "Run an analysis with caches kept to populate these.";
  renderStages(null);
}

function renderResult(outcome) {
  lastOutcome = outcome;
  renderStats([
    ["Samples", outcome.samples.length],
    ["Features", outcome.summary.features],
    ["Identified", outcome.summary.identified],
    ["Compounds", outcome.summary.compounds],
    ["Library", outcome.libraryEntries.toLocaleString()],
    ["Elapsed", `${outcome.elapsedSeconds.toFixed(1)}s`],
  ]);
  $("result").classList.remove("is-hidden");
  $("run-title").textContent = "Analysis complete";
  $("run-lede").textContent = `Finished in ${outcome.elapsedSeconds.toFixed(
    1
  )} s — ${outcome.polarity} mode, ${outcome.summary.features} consensus features.`;
  $("btn-run").textContent = "Run again";
  $("viz-status").textContent = params.keepCache
    ? "Caches are present for this run."
    : "Caches were discarded; re-run with them kept.";
}

function renderPreviousRun(previous) {
  lastOutcome = { outputDir: previous.outputDir };
  renderStats([
    ["Samples", previous.samples],
    ["Features", previous.summary.features],
    ["Identified", previous.summary.identified],
    ["Compounds", previous.summary.compounds],
    ["Reports", previous.reports.length],
    ["Caches", previous.cachesPresent ? "Ready" : "Not kept"],
  ]);
  $("result").classList.remove("is-hidden");
  $("run-title").textContent = "Already processed";
  $("run-lede").textContent = `A completed ${previous.polarity || "DDA"} run is already present in the results folder.`;
  $("btn-run").textContent = "Run again";
  $("viz-status").textContent = previous.cachesPresent
    ? "Existing scan caches are ready for the visualizer."
    : "This run has no scan caches; run again with caches kept for the visualizer.";
  renderStages("complete");
}

/* ------------------------------------------------------------- run cycle */

function setRunning(on) {
  running = on;
  $("btn-run").classList.toggle("is-hidden", on);
  $("btn-cancel").classList.toggle("is-hidden", !on);
  $$("input, .seg button, .btn:not(#btn-cancel)").forEach((el) => {
    if (el.id === "btn-cancel") return;
    el.disabled = on;
  });
  if (visualizer) visualizer.setRunning(on);
  if (!on) {
    $("intensityCutoff").disabled = !$("intensityCutoffOn").checked;
    revalidate();
  }
}

async function startRun() {
  const problems = await revalidate();
  if (problems.some((p) => p.fatal)) {
    showView("run");
    return;
  }
  $("result").classList.add("is-hidden");
  $("run-title").textContent = "Running analysis";
  $("run-lede").textContent = "Running…";
  renderStages("preparing");
  showView("run");
  setRunning(true);
  try {
    await invoke("start_run", { params: collect() });
  } catch (err) {
    setRunning(false);
    log(String(err), "error");
  }
}

function wireEvents() {
  if (!listen) return;

  listen("mk://event", ({ payload }) => {
    switch (payload.type) {
      case "stage":
        renderStages(payload.stage);
        $("rail-stage").textContent = payload.label;
        log(payload.label, "done");
        break;
      case "sample":
        $("rail-stage").textContent = `${payload.index + 1}/${payload.total} ${
          payload.name
        }`;
        break;
      case "progress":
        if (payload.total > 0) {
          $("rail-meter").style.width = `${Math.round(
            (payload.done / payload.total) * 100
          )}%`;
        }
        break;
      case "log":
        log(payload.message, payload.level);
        break;
      case "metric":
        log(`${payload.key}: ${payload.value}`);
        break;
      default:
        break;
    }
  });

  listen("mk://finished", ({ payload }) => {
    setRunning(false);
    $("rail-stage").textContent = "Finished";
    $("rail-meter").style.width = "100%";
    renderStages("complete");
    renderResult(payload);
    if (visualizer) {
      visualizer
        .setOutput(payload.outputDir, { reload: true })
        .catch((err) => log(String(err), "error"));
    }
    log("run complete", "done");
  });

  listen("mk://cancelled", () => {
    setRunning(false);
    $("rail-stage").textContent = "Stopped";
    $("rail-meter").style.width = "0%";
    $("run-lede").textContent = "Stopped before finishing.";
    $("run-title").textContent = "Execution stopped";
    log("run stopped", "warn");
  });

  listen("mk://failed", ({ payload }) => {
    setRunning(false);
    $("rail-stage").textContent = "Failed";
    $("rail-meter").style.width = "0%";
    $("run-lede").textContent = String(payload);
    $("run-title").textContent = "Run failed";
    log(String(payload), "error");
  });
}

/* -------------------------------------------------------------- wiring   */

function wireControls() {
  $$("[data-view]").forEach((btn) => {
    btn.addEventListener("click", () => showView(btn.dataset.view));
  });

  $$(".seg").forEach((seg) => {
    seg.addEventListener("click", (e) => {
      const btn = e.target.closest("button");
      if (!btn || btn.disabled) return;
      Array.from(seg.children).forEach((b) =>
        b.classList.toggle("is-on", b === btn)
      );
      const key = seg.dataset.param;
      if (key && params) {
        params[key] = btn.dataset.value;
        revalidate();
      }
    });
  });

  document.addEventListener("change", (e) => {
    if (e.target.matches("input")) {
      if (e.target.id === "intensityCutoffOn") {
        $("intensityCutoff").disabled = !e.target.checked;
      }
      revalidate();
    }
  });

  $("btn-choose").addEventListener("click", chooseDataset);
  $("btn-rechoose").addEventListener("click", chooseDataset);
  $("btn-rescan").addEventListener("click", () => {
    if (scan) runScan(scan.root);
  });

  $("btn-add-msp").addEventListener("click", () => addLibraries("msp"));
  $("btn-add-csv").addEventListener("click", () => addLibraries("csv"));
  $("btn-libs-dir").addEventListener("click", chooseLibsDir);
  $("btn-libs-dir-2").addEventListener("click", chooseLibsDir);

  $("btn-output-dir").addEventListener("click", async () => {
    const dir = asPath(
      await dialogOpen({ directory: true, title: "Output directory" })
    );
    if (!dir) return;
    $("outputDir").value = dir;
    if (scan) {
      scan.outputDir = dir;
      $("ds-output").textContent = dir;
    }
    if (visualizer) visualizer.setOutput(dir).catch(() => {});
    revalidate();
  });

  $("btn-load-preset").addEventListener("click", async () => {
    const path = asPath(
      await dialogOpen({
        title: "Load settings",
        filters: [{ name: "Settings", extensions: ["json", "txt"] }],
      })
    );
    if (!path) return;
    try {
      params = await invoke("load_preset", { path });
      render();
      await revalidate();
      log(`loaded ${baseName(path)}`);
    } catch (err) {
      log(String(err), "error");
    }
  });

  $("btn-save-preset").addEventListener("click", async () => {
    const path = asPath(
      await dialogSave({
        title: "Save settings",
        defaultPath: "metabokit-settings.json",
        filters: [{ name: "Settings", extensions: ["json"] }],
      })
    );
    if (!path) return;
    try {
      await invoke("save_preset", { path, params: collect() });
      log(`saved ${baseName(path)}`);
    } catch (err) {
      log(String(err), "error");
    }
  });

  $("btn-reset").addEventListener("click", async () => {
    // Reset the analysis settings but keep the dataset the scan established.
    const fresh = await invoke("default_params");
    if (params) {
      fresh.mzmlFiles = params.mzmlFiles;
      fresh.outputDir = params.outputDir;
      fresh.libsDir = params.libsDir;
      fresh.polarity = params.polarity;
      fresh.libraries = params.libraries;
    }
    params = fresh;
    render();
    await revalidate();
    log("analysis settings reset to defaults");
  });

  $("btn-run").addEventListener("click", startRun);
  $("btn-cancel").addEventListener("click", () => {
    invoke("cancel_run").catch(() => {});
    $("rail-stage").textContent = "Stopping…";
  });
  $("btn-clear-log").addEventListener("click", () => clear($("log")));
  $("btn-open-output").addEventListener("click", async () => {
    const path = lastOutcome ? lastOutcome.outputDir : params && params.outputDir;
    await openOutputDirectory(path);
  });

  $("theme-seg").addEventListener("click", (e) => {
    const btn = e.target.closest("button");
    if (!btn) return;
    document.documentElement.dataset.theme = btn.dataset.value;
    try {
      localStorage.setItem("mk-theme", btn.dataset.value);
    } catch (err) {
      /* private mode; the choice simply will not persist */
    }
  });
}

/* --------------------------------------------------------------- startup */

async function main() {
  let stored = "system";
  try {
    stored = localStorage.getItem("mk-theme") || "system";
  } catch (err) {
    /* ignore */
  }
  document.documentElement.dataset.theme = stored;
  Array.from($("theme-seg").children).forEach((b) =>
    b.classList.toggle("is-on", b.dataset.value === stored)
  );

  renderStages(null);
  visualizer = createVisualizer({ invoke, log });
  wireControls();
  wireEvents();

  if (!invoke) {
    log("running without the Tauri bridge — controls are inert", "warn");
    return;
  }

  params = await invoke("default_params");
  render();
  await refreshSystem();
  await revalidate();
  if (await invoke("is_running")) setRunning(true);

  const previousDataset = storedDataset();
  if (previousDataset) {
    let stillExists = false;
    try {
      stillExists = await invoke("directory_exists", { path: previousDataset });
    } catch (err) {
      log(`could not check the previous dataset: ${err}`, "warn");
    }

    if (stillExists) {
      await runScan(previousDataset, { restored: true });
    } else {
      rememberDataset(null);
      $("rail-stage").textContent = "Choose a dataset";
      log(`previous dataset is no longer available: ${previousDataset}`, "warn");
    }
  }
}

main();
