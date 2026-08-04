/* Integrated MetaboKit visualizer.
 *
 * Canvas is used instead of one SVG node per point/peak. The feature map has
 * a typed-array spatial grid for O(nearby points) hit testing, redraws are
 * requestAnimationFrame-coalesced, and mirror canvases are materialized only
 * while near the viewport.
 */

const $ = (id) => document.getElementById(id);
const IS_MAC = /Macintosh|Mac OS X/.test(navigator.userAgent);
if (IS_MAC) document.documentElement.classList.add("is-mac");
const fmt = (value, digits = 2) =>
  Number.isFinite(value) ? Number(value).toFixed(digits) : "—";

function bytes(value) {
  if (!value) return "0 B";
  const units = ["B", "KB", "MB", "GB"];
  let n = value;
  let at = 0;
  while (n >= 1024 && at < units.length - 1) {
    n /= 1024;
    at += 1;
  }
  return `${n < 10 && at ? n.toFixed(1) : Math.round(n)} ${units[at]}`;
}

function palette() {
  const style = getComputedStyle(document.body);
  return {
    paper: style.getPropertyValue("--paper").trim(),
    sunk: style.getPropertyValue("--paper-sunk").trim(),
    ink: style.getPropertyValue("--ink").trim(),
    ink2: style.getPropertyValue("--ink-2").trim(),
    ink3: style.getPropertyValue("--ink-3").trim(),
    rule: style.getPropertyValue("--rule").trim(),
    accent: style.getPropertyValue("--accent").trim(),
    bad: style.getPropertyValue("--bad").trim(),
    good: style.getPropertyValue("--good").trim(),
    font: style.getPropertyValue("--sans").trim(),
    mono: style.getPropertyValue("--mono").trim(),
  };
}

function prepareCanvas(canvas, dprCap = 2) {
  const rect = canvas.getBoundingClientRect();
  const width = Math.max(1, Math.round(rect.width));
  const height = Math.max(1, Math.round(rect.height));
  const dpr = Math.min(window.devicePixelRatio || 1, dprCap);
  const pixelWidth = Math.max(1, Math.round(width * dpr));
  const pixelHeight = Math.max(1, Math.round(height * dpr));
  if (canvas.width !== pixelWidth || canvas.height !== pixelHeight) {
    canvas.width = pixelWidth;
    canvas.height = pixelHeight;
  }
  const ctx = canvas.getContext("2d", { alpha: false, desynchronized: true });
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, width, height };
}

function ticks(lo, hi, count = 5) {
  if (!Number.isFinite(lo) || !Number.isFinite(hi) || hi <= lo) return [lo || 0];
  const rough = (hi - lo) / Math.max(1, count);
  const power = 10 ** Math.floor(Math.log10(rough));
  const ratio = rough / power;
  const step = (ratio >= 5 ? 10 : ratio >= 2 ? 5 : ratio >= 1 ? 2 : 1) * power;
  const first = Math.ceil(lo / step) * step;
  const out = [];
  for (let value = first; value <= hi + step * 0.01; value += step) out.push(value);
  return out;
}

function plotAxes(ctx, box, xDomain, yDomain, labels = {}) {
  const p = palette();
  const x = (value) =>
    box.left + ((value - xDomain[0]) / (xDomain[1] - xDomain[0] || 1)) * box.width;
  const y = (value) =>
    box.top + box.height - ((value - yDomain[0]) / (yDomain[1] - yDomain[0] || 1)) * box.height;
  ctx.lineWidth = 1;
  ctx.font = `11px ${p.mono}`;
  ctx.fillStyle = p.ink3;
  ctx.strokeStyle = p.rule;
  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  for (const value of ticks(xDomain[0], xDomain[1], 6)) {
    const at = x(value);
    ctx.beginPath();
    ctx.moveTo(at, box.top);
    ctx.lineTo(at, box.top + box.height);
    ctx.stroke();
    ctx.fillText(formatTick(value), at, box.top + box.height + 8);
  }
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (const value of ticks(yDomain[0], yDomain[1], 5)) {
    const at = y(value);
    ctx.beginPath();
    ctx.moveTo(box.left, at);
    ctx.lineTo(box.left + box.width, at);
    ctx.stroke();
    ctx.fillText(formatTick(value), box.left - 9, at);
  }
  ctx.fillStyle = p.ink2;
  if (labels.x) {
    ctx.textAlign = "right";
    ctx.textBaseline = "bottom";
    ctx.fillText(labels.x, box.left + box.width, box.top + box.height + 37);
  }
  if (labels.y) {
    ctx.save();
    ctx.translate(13, box.top);
    ctx.textAlign = "left";
    ctx.textBaseline = "top";
    ctx.fillText(labels.y, 0, 0);
    ctx.restore();
  }
  return { x, y };
}

function formatTick(value) {
  const abs = Math.abs(value);
  if ((abs >= 10_000 || (abs > 0 && abs < 0.001))) return value.toExponential(1);
  if (abs >= 100) return value.toFixed(0);
  if (abs >= 10) return value.toFixed(1);
  return value.toFixed(2).replace(/0+$/, "").replace(/\.$/, "");
}

function showTip(element, wrap, x, y, text) {
  element.textContent = text;
  element.classList.remove("is-hidden");
  const maxX = Math.max(8, wrap.clientWidth - element.offsetWidth - 8);
  const maxY = Math.max(8, wrap.clientHeight - element.offsetHeight - 8);
  element.style.left = `${Math.max(8, Math.min(maxX, x + 13))}px`;
  element.style.top = `${Math.max(8, Math.min(maxY, y - element.offsetHeight - 10))}px`;
}

function hideTip(element) {
  element.classList.add("is-hidden");
}

function extent(values) {
  let lo = Infinity;
  let hi = -Infinity;
  for (const value of values) {
    if (!Number.isFinite(value)) continue;
    if (value < lo) lo = value;
    if (value > hi) hi = value;
  }
  if (!Number.isFinite(lo)) return [0, 1];
  if (lo === hi) return [lo - 0.5, hi + 0.5];
  return [lo, hi];
}

function padded([lo, hi], fraction = 0.025) {
  const pad = (hi - lo || 1) * fraction;
  return [lo - pad, hi + pad];
}

function compoundKey(name, fallback) {
  const cleaned = String(name || "")
    .replace(/;\s*CE\s*[-+]?\d+(?:\.\d+)?(?=;|$)/gi, "")
    .replace(/\s*;\s*/g, "; ")
    .trim();
  return cleaned || fallback;
}

export function createVisualizer({ invoke, log }) {
  const els = {
    sample: $("viz-sample"),
    reload: $("viz-reload"),
    reset: $("viz-reset"),
    intensityMode: $("viz-intensity-mode"),
    status: $("viz-status"),
    empty: $("viz-empty"),
    emptyTitle: $("viz-empty-title"),
    emptyDetail: $("viz-empty-detail"),
    workspace: $("viz-workspace"),
    facts: $("viz-facts"),
    map: $("viz-map"),
    mapTip: $("viz-map-tip"),
    selection: $("viz-selection"),
    featureTitle: $("viz-feature-title"),
    featureCoords: $("viz-feature-coords"),
    xic: $("viz-xic"),
    xicTip: $("viz-xic-tip"),
    xicMeta: $("viz-xic-meta"),
    spectrum: $("viz-spectrum"),
    spectrumTip: $("viz-spectrum-tip"),
    spectrumTitle: $("viz-spectrum-title"),
    spectrumMeta: $("viz-spectrum-meta"),
    mirrorCount: $("viz-mirror-count"),
    mirrorList: $("viz-mirror-list"),
  };

  let outputDir = null;
  let active = false;
  let analysisRunning = false;
  let session = null;
  let overview = null;
  let selected = null;
  let detail = null;
  let currentSpectrum = null;
  let loadToken = 0;
  let mapBase = null;
  let mapDomain = null;
  let mapGrid = null;
  let mapFrame = 0;
  let spectrumBase = null;
  let spectrumDomain = null;
  let mirrorObserver = null;
  let intensityMode = "raw";
  try {
    intensityMode = localStorage.getItem("mk-viz-intensity") === "sqrt" ? "sqrt" : "raw";
  } catch (error) {
    /* storage is optional */
  }
  Array.from(els.intensityMode.children).forEach((button) =>
    button.classList.toggle("is-on", button.dataset.value === intensityMode)
  );

  const transformedIntensity = (value) =>
    intensityMode === "sqrt" ? Math.sqrt(Math.max(0, value)) : Math.max(0, value);
  const detailDpr = () => (IS_MAC ? 1.25 : 1.5);
  const mirrorDpr = () => (IS_MAC ? 1.15 : 1.35);

  function status(text) {
    els.status.textContent = text;
  }

  function empty(title, detailText) {
    els.emptyTitle.textContent = title;
    els.emptyDetail.textContent = detailText;
    els.empty.classList.remove("is-hidden");
    els.workspace.classList.add("is-hidden");
    els.selection.classList.add("is-hidden");
  }

  function clearMirrors() {
    if (mirrorObserver) mirrorObserver.disconnect();
    mirrorObserver = null;
    els.mirrorList.textContent = "";
    els.mirrorCount.textContent = "";
  }

  async function loadSession(force = false) {
    if (!active || !invoke || !outputDir || analysisRunning) return;
    const token = ++loadToken;
    status("Reading cached run…");
    els.reload.disabled = true;
    try {
      const next = await invoke("visualizer_open", { outputDir });
      if (token !== loadToken) return;
      session = next;
      const previous = force ? null : els.sample.value;
      els.sample.textContent = "";
      next.samples.forEach((sample) => {
        const option = document.createElement("option");
        option.value = sample.id;
        option.textContent = sample.label;
        option.title = `${sample.label} · ${bytes(sample.cacheBytes)}`;
        els.sample.append(option);
      });
      if (!next.samples.length) {
        els.sample.disabled = true;
        els.reset.disabled = true;
        empty(
          "No scan caches found",
          "Run this dataset again with Keep scan caches enabled. Reports alone can restore the run summary, but raw chromatograms and spectra require the caches.",
        );
        status("No visualizer caches in this results folder.");
        return;
      }
      if (previous && next.samples.some((sample) => sample.id === previous)) {
        els.sample.value = previous;
      }
      els.sample.disabled = false;
      els.reload.disabled = false;
      els.reset.disabled = false;
      await loadOverview(token);
    } catch (error) {
      if (token !== loadToken) return;
      empty("Could not open visualizer data", String(error));
      status("Visualizer data could not be opened.");
      log(String(error), "error");
      els.reload.disabled = false;
    }
  }

  async function loadOverview(parentToken = ++loadToken) {
    if (!outputDir || !els.sample.value || analysisRunning) return;
    const token = parentToken;
    const sample = els.sample.value;
    status(`Loading ${sample}…`);
    els.sample.disabled = true;
    els.selection.classList.add("is-hidden");
    clearMirrors();
    selected = null;
    detail = null;
    try {
      const next = await invoke("visualizer_overview", { outputDir, sample });
      if (token !== loadToken) return;
      overview = next;
      const rt = padded(extent(next.rt));
      const mz = padded(extent(next.mz));
      mapBase = { rt, mz };
      mapDomain = { rt: [...rt], mz: [...mz] };
      renderFacts();
      els.empty.classList.add("is-hidden");
      els.workspace.classList.remove("is-hidden");
      els.sample.disabled = false;
      const shown = next.mz.length.toLocaleString();
      const total = next.total.toLocaleString();
      status(
        next.truncated
          ? `${shown} of ${total} features shown · bounded overview sample`
          : `${total} features · ${next.source}`,
      );
      requestMapDraw();
    } catch (error) {
      if (token !== loadToken) return;
      empty("Could not load the feature map", String(error));
      status("Feature map failed to load.");
      log(String(error), "error");
      els.sample.disabled = false;
    }
  }

  function renderFacts() {
    els.facts.textContent = "";
    const identified = overview.names.reduce((sum, name) => sum + (name ? 1 : 0), 0);
    const selectedSample = session.samples.find((item) => item.id === els.sample.value);
    const facts = [
      ["Plotted", overview.mz.length.toLocaleString()],
      ["Identified", identified.toLocaleString()],
      ["Cache", bytes(selectedSample ? selectedSample.cacheBytes : 0)],
      ["Source", overview.source],
    ];
    for (const [label, value] of facts) {
      const wrap = document.createElement("div");
      const dt = document.createElement("dt");
      const dd = document.createElement("dd");
      dt.textContent = label;
      dd.textContent = value;
      dd.title = value;
      wrap.append(dt, dd);
      els.facts.append(wrap);
    }
  }

  function requestMapDraw() {
    if (mapFrame) return;
    mapFrame = requestAnimationFrame(() => {
      mapFrame = 0;
      drawMap();
    });
  }

  function drawMap() {
    if (!overview || !mapDomain || els.workspace.classList.contains("is-hidden")) return;
    const { ctx, width, height } = prepareCanvas(els.map, IS_MAC ? 1.5 : 2);
    if (width < 40 || height < 40) return;
    const p = palette();
    ctx.fillStyle = p.paper;
    ctx.fillRect(0, 0, width, height);
    const box = { left: 62, top: 24, width: Math.max(1, width - 84), height: Math.max(1, height - 76) };
    const scale = plotAxes(ctx, box, mapDomain.rt, mapDomain.mz, { x: "retention time (min)", y: "m/z" });

    const cell = 18;
    const cols = Math.ceil(box.width / cell);
    const rows = Math.ceil(box.height / cell);
    const heads = new Int32Array(Math.max(1, cols * rows));
    heads.fill(-1);
    const next = new Int32Array(overview.mz.length);
    next.fill(-1);
    const identified = [];
    const unknown = [];
    for (let i = 0; i < overview.mz.length; i += 1) {
      const x = scale.x(overview.rt[i]);
      const y = scale.y(overview.mz[i]);
      if (x < box.left || x > box.left + box.width || y < box.top || y > box.top + box.height) continue;
      (overview.names[i] ? identified : unknown).push([x, y]);
      const col = Math.max(0, Math.min(cols - 1, Math.floor((x - box.left) / cell)));
      const row = Math.max(0, Math.min(rows - 1, Math.floor((y - box.top) / cell)));
      const at = row * cols + col;
      next[i] = heads[at];
      heads[at] = i;
    }
    ctx.lineWidth = 1;
    ctx.strokeStyle = p.ink2;
    ctx.beginPath();
    for (const [x, y] of unknown) {
      ctx.moveTo(x + 2.4, y);
      ctx.arc(x, y, 2.4, 0, Math.PI * 2);
    }
    ctx.stroke();
    ctx.fillStyle = p.accent;
    ctx.beginPath();
    for (const [x, y] of identified) {
      ctx.moveTo(x + 3, y);
      ctx.arc(x, y, 3, 0, Math.PI * 2);
    }
    ctx.fill();

    if (selected) {
      const x = scale.x(selected.rt);
      const y = scale.y(selected.mz);
      const pointRadius = selected.name ? 3 : 2.4;
      ctx.strokeStyle = p.ink;
      ctx.lineWidth = 2;
      ctx.beginPath();
      // With a 2 px stroke, radius + 1 puts the inner edge exactly against
      // the selected point rather than leaving a halo of empty space.
      ctx.arc(x, y, pointRadius + 1, 0, Math.PI * 2);
      ctx.stroke();
    }
    mapGrid = { box, scale, cell, cols, rows, heads, next };
  }

  function nearestMapPoint(x, y, radius = 15) {
    if (!mapGrid) return -1;
    const { box, cell, cols, rows, heads, next, scale } = mapGrid;
    if (x < box.left || x > box.left + box.width || y < box.top || y > box.top + box.height) return -1;
    const col = Math.floor((x - box.left) / cell);
    const row = Math.floor((y - box.top) / cell);
    let best = -1;
    let bestDistance = radius * radius;
    for (let yy = Math.max(0, row - 1); yy <= Math.min(rows - 1, row + 1); yy += 1) {
      for (let xx = Math.max(0, col - 1); xx <= Math.min(cols - 1, col + 1); xx += 1) {
        for (let i = heads[yy * cols + xx]; i >= 0; i = next[i]) {
          const dx = scale.x(overview.rt[i]) - x;
          const dy = scale.y(overview.mz[i]) - y;
          const distance = dx * dx + dy * dy;
          if (distance < bestDistance) {
            best = i;
            bestDistance = distance;
          }
        }
      }
    }
    return best;
  }

  function mapCoordinates(event) {
    const rect = els.map.getBoundingClientRect();
    return [event.clientX - rect.left, event.clientY - rect.top];
  }

  let mapDrag = null;
  els.map.addEventListener("pointerdown", (event) => {
    if (!mapDomain || !mapGrid) return;
    const [x, y] = mapCoordinates(event);
    mapDrag = { x, y, moved: false, rt: [...mapDomain.rt], mz: [...mapDomain.mz] };
    els.map.setPointerCapture(event.pointerId);
  });
  els.map.addEventListener("pointermove", (event) => {
    const [x, y] = mapCoordinates(event);
    if (mapDrag && (event.buttons & 1)) {
      const dx = x - mapDrag.x;
      const dy = y - mapDrag.y;
      if (Math.abs(dx) + Math.abs(dy) > 3) mapDrag.moved = true;
      const rtSpan = mapDrag.rt[1] - mapDrag.rt[0];
      const mzSpan = mapDrag.mz[1] - mapDrag.mz[0];
      const rtShift = (-dx / mapGrid.box.width) * rtSpan;
      const mzShift = (dy / mapGrid.box.height) * mzSpan;
      mapDomain.rt = constrainDomain(
        [mapDrag.rt[0] + rtShift, mapDrag.rt[1] + rtShift],
        mapBase.rt,
      );
      mapDomain.mz = constrainDomain(
        [mapDrag.mz[0] + mzShift, mapDrag.mz[1] + mzShift],
        mapBase.mz,
      );
      hideTip(els.mapTip);
      requestMapDraw();
      return;
    }
    const i = nearestMapPoint(x, y);
    if (i < 0) {
      hideTip(els.mapTip);
      return;
    }
    const name = overview.names[i] || "Unidentified feature";
    const quality = overview.source === "feature cache" ? "smoothness" : "S/N";
    const metrics = `m/z ${fmt(overview.mz[i], 4)} · RT ${fmt(overview.rt[i], 3)} min\nshape ${fmt(overview.shape[i])} · ${quality} ${fmt(overview.smoothness[i])}`;
    showTip(els.mapTip, els.map.parentElement, x, y, `${name}\n${metrics}`);
  });
  els.map.addEventListener("pointerup", (event) => {
    if (!mapDrag) return;
    const [x, y] = mapCoordinates(event);
    const moved = mapDrag.moved;
    mapDrag = null;
    if (!moved) {
      const i = nearestMapPoint(x, y);
      if (i >= 0) selectFeature(i);
    }
  });
  els.map.addEventListener("pointerleave", () => hideTip(els.mapTip));
  els.map.addEventListener("wheel", (event) => {
    if (!mapGrid || !mapDomain) return;
    event.preventDefault();
    const [x, y] = mapCoordinates(event);
    const factor = Math.exp(Math.max(-1, Math.min(1, event.deltaY * 0.0015)));
    const centerRt = mapDomain.rt[0] + ((x - mapGrid.box.left) / mapGrid.box.width) * (mapDomain.rt[1] - mapDomain.rt[0]);
    const centerMz = mapDomain.mz[1] - ((y - mapGrid.box.top) / mapGrid.box.height) * (mapDomain.mz[1] - mapDomain.mz[0]);
    mapDomain.rt = zoomDomain(mapDomain.rt, centerRt, factor, mapBase.rt);
    mapDomain.mz = zoomDomain(mapDomain.mz, centerMz, factor, mapBase.mz);
    hideTip(els.mapTip);
    requestMapDraw();
  }, { passive: false });
  els.map.addEventListener("dblclick", resetMap);

  function zoomDomain(domain, center, factor, base) {
    const minSpan = (base[1] - base[0]) / 999;
    let lo = center + (domain[0] - center) * factor;
    let hi = center + (domain[1] - center) * factor;
    if (hi - lo < minSpan) {
      lo = center - minSpan / 2;
      hi = center + minSpan / 2;
    }
    if (hi - lo >= base[1] - base[0]) return [...base];
    return constrainDomain([lo, hi], base);
  }

  function constrainDomain(domain, base) {
    const span = domain[1] - domain[0];
    const baseSpan = base[1] - base[0];
    if (span >= baseSpan) return [...base];
    let [lo, hi] = domain;
    if (lo < base[0]) {
      hi += base[0] - lo;
      lo = base[0];
    }
    if (hi > base[1]) {
      lo -= hi - base[1];
      hi = base[1];
    }
    return [lo, hi];
  }

  function resetMap() {
    if (!mapBase) return;
    mapDomain = { rt: [...mapBase.rt], mz: [...mapBase.mz] };
    hideTip(els.mapTip);
    requestMapDraw();
  }

  async function selectFeature(index) {
    if (analysisRunning) return;
    const token = ++loadToken;
    selected = {
      index,
      mz: overview.mz[index],
      rt: overview.rt[index],
      halfWidth: overview.halfWidth[index],
      shape: overview.shape[index],
      smoothness: overview.smoothness[index],
      name: overview.names[index],
      group: overview.group[index],
    };
    requestMapDraw();
    els.selection.classList.remove("is-hidden");
    els.featureTitle.textContent = selected.name || "Unidentified feature";
    els.featureCoords.textContent = `${fmt(selected.mz, 4)} m/z · ${fmt(selected.rt, 3)} min`;
    els.xicMeta.textContent = "Loading…";
    els.spectrumMeta.textContent = "";
    clearMirrors();
    drawLoading(els.xic, "Extracting chromatogram…");
    drawLoading(els.spectrum, "Loading spectra…");
    try {
      const next = await invoke("visualizer_feature", {
        outputDir,
        sample: els.sample.value,
        mz: selected.mz,
        rt: selected.rt,
        halfWidth: selected.halfWidth,
      });
      if (token !== loadToken) return;
      detail = next;
      detail.fragmentAt = assignFragmentation(next.chromatogram, next.spectra);
      els.xicMeta.textContent = `${next.chromatogram.length} MS1 scans · ${next.spectra.length} MS/MS events${next.spectraTruncated ? " (nearest shown)" : ""}`;
      currentSpectrum = nearestSpectrum(next.spectra, selected.rt);
      spectrumBase = currentSpectrum ? spectrumExtent(currentSpectrum) : [0, 1];
      spectrumDomain = [...spectrumBase];
      drawXic();
      drawSpectrum();
      renderMirrors(next.mirrors);
    } catch (error) {
      if (token !== loadToken) return;
      els.xicMeta.textContent = "Could not load this feature";
      drawLoading(els.xic, String(error));
      drawLoading(els.spectrum, "No spectrum available");
      log(String(error), "error");
    }
  }

  function drawLoading(canvas, message) {
    const { ctx, width, height } = prepareCanvas(canvas, detailDpr());
    const p = palette();
    ctx.fillStyle = p.paper;
    ctx.fillRect(0, 0, width, height);
    ctx.fillStyle = p.ink3;
    ctx.font = `12px ${p.mono}`;
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(message, width / 2, height / 2);
  }

  function nearestSpectrum(spectra, rt) {
    if (!spectra || !spectra.length) return null;
    return spectra.reduce((best, item) =>
      Math.abs(item.rt - rt) < Math.abs(best.rt - rt) ? item : best,
    );
  }

  function assignFragmentation(chromatogram, spectra) {
    const assigned = new Int32Array(chromatogram.length);
    assigned.fill(-1);
    if (!chromatogram.length) return assigned;
    for (let spectrumIndex = 0; spectrumIndex < spectra.length; spectrumIndex += 1) {
      const rt = spectra[spectrumIndex].rt;
      let lo = 0;
      let hi = chromatogram.length;
      while (lo < hi) {
        const middle = (lo + hi) >>> 1;
        if (chromatogram[middle][0] < rt) lo = middle + 1;
        else hi = middle;
      }
      const right = Math.min(chromatogram.length - 1, lo);
      const left = Math.max(0, right - 1);
      const index = Math.abs(chromatogram[left][0] - rt) <= Math.abs(chromatogram[right][0] - rt)
        ? left
        : right;
      const previous = assigned[index];
      if (previous < 0 || Math.abs(rt - chromatogram[index][0]) < Math.abs(spectra[previous].rt - chromatogram[index][0])) {
        assigned[index] = spectrumIndex;
      }
    }
    return assigned;
  }

  function drawXic() {
    if (!detail || !selected) return;
    const { ctx, width, height } = prepareCanvas(els.xic, detailDpr());
    const p = palette();
    ctx.fillStyle = p.paper;
    ctx.fillRect(0, 0, width, height);
    if (!detail.chromatogram.length) {
      drawLoading(els.xic, "No MS1 points in this window");
      return;
    }
    const xDomain = extent(detail.chromatogram.map((point) => point[0]));
    const yDomain = [0, Math.max(1, ...detail.chromatogram.map((point) => point[1])) * 1.06];
    const box = { left: 60, top: 19, width: Math.max(1, width - 78), height: Math.max(1, height - 58) };
    const scale = plotAxes(ctx, box, xDomain, yDomain, { x: "RT (min)", y: "intensity" });
    const half = Math.max(0, detail.integrationHalfWidth || 0);
    ctx.fillStyle = `${p.accent}18`;
    ctx.fillRect(scale.x(selected.rt - half), box.top, Math.max(1, scale.x(selected.rt + half) - scale.x(selected.rt - half)), box.height);
    ctx.strokeStyle = p.ink;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    detail.chromatogram.forEach((point, i) => {
      const x = scale.x(point[0]);
      const y = scale.y(point[1]);
      if (i) ctx.lineTo(x, y);
      else ctx.moveTo(x, y);
    });
    ctx.stroke();

    // Every chromatographic scan remains visible as a point. Fragmented scans
    // are assigned to their nearest MS1 acquisition and drawn as a larger,
    // coloured point so scan density and DDA coverage are readable together.
    const fragmentAt = detail.fragmentAt || assignFragmentation(detail.chromatogram, detail.spectra);
    ctx.fillStyle = p.ink3;
    ctx.beginPath();
    detail.chromatogram.forEach((point, index) => {
      if (fragmentAt[index] >= 0) return;
      ctx.moveTo(scale.x(point[0]) + 1.8, scale.y(point[1]));
      ctx.arc(scale.x(point[0]), scale.y(point[1]), 1.8, 0, Math.PI * 2);
    });
    ctx.fill();
    ctx.fillStyle = p.bad;
    ctx.beginPath();
    detail.chromatogram.forEach((point, index) => {
      if (fragmentAt[index] < 0) return;
      ctx.moveTo(scale.x(point[0]) + 3.2, scale.y(point[1]));
      ctx.arc(scale.x(point[0]), scale.y(point[1]), 3.2, 0, Math.PI * 2);
    });
    ctx.fill();

    ctx.strokeStyle = `${p.bad}88`;
    ctx.lineWidth = 1;
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    for (const spectrum of detail.spectra) {
      if (currentSpectrum && spectrum.id === currentSpectrum.id) continue;
      ctx.moveTo(scale.x(spectrum.rt), box.top);
      ctx.lineTo(scale.x(spectrum.rt), box.top + box.height);
    }
    ctx.stroke();
    if (currentSpectrum) {
      ctx.strokeStyle = p.bad;
      ctx.lineWidth = 2.5;
      ctx.beginPath();
      ctx.moveTo(scale.x(currentSpectrum.rt), box.top);
      ctx.lineTo(scale.x(currentSpectrum.rt), box.top + box.height);
      ctx.stroke();
    }
    ctx.setLineDash([]);
    els.xic._viz = { box, scale, xDomain };
  }

  els.xic.addEventListener("pointermove", (event) => {
    if (!detail || !detail.chromatogram.length || !els.xic._viz) return;
    const rect = els.xic.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;
    const { box, xDomain } = els.xic._viz;
    const rt = xDomain[0] + ((x - box.left) / box.width) * (xDomain[1] - xDomain[0]);
    let lo = 0;
    let hi = detail.chromatogram.length;
    while (lo < hi) {
      const middle = (lo + hi) >>> 1;
      if (detail.chromatogram[middle][0] < rt) lo = middle + 1;
      else hi = middle;
    }
    const right = Math.min(detail.chromatogram.length - 1, lo);
    const left = Math.max(0, right - 1);
    const index = Math.abs(detail.chromatogram[left][0] - rt) <= Math.abs(detail.chromatogram[right][0] - rt)
      ? left
      : right;
    const point = detail.chromatogram[index];
    const spectrumIndex = detail.fragmentAt ? detail.fragmentAt[index] : -1;
    const spectrum = spectrumIndex >= 0 ? detail.spectra[spectrumIndex] : null;
    if (spectrum && (!currentSpectrum || spectrum.id !== currentSpectrum.id)) {
      currentSpectrum = spectrum;
      spectrumBase = spectrumExtent(currentSpectrum);
      spectrumDomain = [...spectrumBase];
      drawXic();
      drawSpectrum();
    }
    showTip(
      els.xicTip,
      els.xic.parentElement,
      x,
      y,
      `MS1 scan ${index + 1}\nRT ${fmt(point[0], 3)} min · ${formatTick(point[1])} intensity${spectrum ? `\nMS/MS acquired · ${fmt(spectrum.precursorMz, 4)} m/z · CE ${fmt(spectrum.ce, 1)}` : "\nNo MS/MS fragmentation"}`,
    );
  });
  els.xic.addEventListener("pointerleave", () => hideTip(els.xicTip));

  function spectrumExtent(spectrum) {
    if (!spectrum || !spectrum.peaks.length) return [0, Math.max(1, spectrum ? spectrum.precursorMz : 1)];
    return [0, Math.max(spectrum.precursorMz, spectrum.peaks[spectrum.peaks.length - 1][0]) * 1.03];
  }

  function drawSpectrum() {
    if (!currentSpectrum) {
      els.spectrumTitle.textContent = "Acquired spectrum";
      els.spectrumMeta.textContent = "No nearby MS/MS scan";
      drawLoading(els.spectrum, "No nearby MS/MS spectrum");
      return;
    }
    const { ctx, width, height } = prepareCanvas(els.spectrum, detailDpr());
    const p = palette();
    ctx.fillStyle = p.paper;
    ctx.fillRect(0, 0, width, height);
    const maxI = Math.max(1, ...currentSpectrum.peaks.map((peak) => transformedIntensity(peak[1])));
    const box = { left: 58, top: 19, width: Math.max(1, width - 76), height: Math.max(1, height - 58) };
    const scale = plotAxes(ctx, box, spectrumDomain, [0, maxI * 1.05], {
      x: "fragment m/z",
      y: intensityMode === "sqrt" ? "√ intensity" : "intensity",
    });
    ctx.strokeStyle = p.ink;
    ctx.lineWidth = 1.3;
    ctx.beginPath();
    for (const peak of currentSpectrum.peaks) {
      const x = scale.x(peak[0]);
      if (x < box.left || x > box.left + box.width) continue;
      ctx.moveTo(x, scale.y(0));
      ctx.lineTo(x, scale.y(transformedIntensity(peak[1])));
    }
    ctx.stroke();
    const precursorX = scale.x(currentSpectrum.precursorMz);
    if (precursorX >= box.left && precursorX <= box.left + box.width) {
      ctx.fillStyle = p.bad;
      ctx.beginPath();
      ctx.moveTo(precursorX, box.top + box.height);
      ctx.lineTo(precursorX - 6, box.top + box.height - 10);
      ctx.lineTo(precursorX + 6, box.top + box.height - 10);
      ctx.closePath();
      ctx.fill();
    }
    els.spectrumTitle.textContent = `MS/MS at ${fmt(currentSpectrum.precursorMz, 4)} m/z`;
    els.spectrumMeta.textContent = `${fmt(currentSpectrum.rt, 3)} min · CE ${fmt(currentSpectrum.ce, 1)}`;
    els.spectrum._viz = { box, scale };
  }

  function spectrumCoordinates(event) {
    const rect = els.spectrum.getBoundingClientRect();
    return [event.clientX - rect.left, event.clientY - rect.top];
  }

  let spectrumDrag = null;
  els.spectrum.addEventListener("pointerdown", (event) => {
    if (!spectrumDomain || !els.spectrum._viz) return;
    spectrumDrag = { x: spectrumCoordinates(event)[0], domain: [...spectrumDomain] };
    els.spectrum.setPointerCapture(event.pointerId);
  });
  els.spectrum.addEventListener("pointermove", (event) => {
    if (!currentSpectrum || !els.spectrum._viz) return;
    const [x, y] = spectrumCoordinates(event);
    if (spectrumDrag && (event.buttons & 1)) {
      const shift = (-(x - spectrumDrag.x) / els.spectrum._viz.box.width) * (spectrumDrag.domain[1] - spectrumDrag.domain[0]);
      spectrumDomain = constrainDomain(
        [spectrumDrag.domain[0] + shift, spectrumDrag.domain[1] + shift],
        spectrumBase,
      );
      hideTip(els.spectrumTip);
      drawSpectrum();
      return;
    }
    let best = null;
    let distance = 10;
    for (const peak of currentSpectrum.peaks) {
      const at = els.spectrum._viz.scale.x(peak[0]);
      const next = Math.abs(at - x);
      if (next < distance) {
        best = peak;
        distance = next;
      }
    }
    if (best) showTip(
      els.spectrumTip,
      els.spectrum.parentElement,
      x,
      y,
      `${fmt(best[0], 4)} m/z\n${formatTick(transformedIntensity(best[1]))} ${intensityMode === "sqrt" ? "√ intensity" : "intensity"}`,
    );
    else hideTip(els.spectrumTip);
  });
  els.spectrum.addEventListener("pointerup", () => { spectrumDrag = null; });
  els.spectrum.addEventListener("pointerleave", () => {
    spectrumDrag = null;
    hideTip(els.spectrumTip);
  });
  els.spectrum.addEventListener("wheel", (event) => {
    if (!spectrumDomain || !spectrumBase || !els.spectrum._viz) return;
    event.preventDefault();
    const x = spectrumCoordinates(event)[0];
    const center = spectrumDomain[0] + ((x - els.spectrum._viz.box.left) / els.spectrum._viz.box.width) * (spectrumDomain[1] - spectrumDomain[0]);
    const factor = Math.exp(Math.max(-1, Math.min(1, event.deltaY * 0.0015)));
    spectrumDomain = zoomDomain(spectrumDomain, center, factor, spectrumBase);
    drawSpectrum();
  }, { passive: false });
  els.spectrum.addEventListener("dblclick", () => {
    if (!spectrumBase) return;
    spectrumDomain = [...spectrumBase];
    drawSpectrum();
  });

  function renderMirrors(mirrors) {
    clearMirrors();
    if (!mirrors.length) {
      els.mirrorCount.textContent = "No library matches";
      const emptyMessage = document.createElement("p");
      emptyMessage.className = "viz-mirror-empty";
      emptyMessage.textContent = "No ranked library spectrum was recorded for this feature in this sample.";
      els.mirrorList.append(emptyMessage);
      return;
    }

    const groups = [];
    const byName = new Map();
    mirrors.forEach((mirror, index) => {
      const key = compoundKey(mirror.name, `Unlabelled match ${index + 1}`);
      let group = byName.get(key);
      if (!group) {
        group = { name: key, mirrors: [] };
        byName.set(key, group);
        groups.push(group);
      }
      group.mirrors.push(mirror);
    });
    els.mirrorCount.textContent = `${mirrors.length} librar${mirrors.length === 1 ? "y spectrum" : "y spectra"} · ${groups.length} compound${groups.length === 1 ? "" : "s"} · best score first`;

    const root = document.querySelector(".stage");
    mirrorObserver = new IntersectionObserver(
      (entries) => {
        for (const entry of entries) {
          const canvas = entry.target;
          if (entry.isIntersecting) {
            canvas.dataset.live = "true";
            drawMirror(canvas, canvas._mirror);
          } else {
            canvas.dataset.live = "false";
            canvas.width = 1;
            canvas.height = 1;
          }
        }
      },
      { root, rootMargin: `${IS_MAC ? 300 : 650}px 0px` },
    );

    groups.forEach((group) => {
      const mirror = group.mirrors[0];
      const card = document.createElement("article");
      card.className = "viz-mirror-card";
      const head = document.createElement("header");
      head.className = "viz-mirror-head";
      const title = document.createElement("h4");
      title.textContent = group.name;
      title.title = group.name;
      const meta = document.createElement("div");
      meta.className = "viz-mirror-meta mono micro";
      const score = document.createElement("span");
      const mass = document.createElement("span");
      meta.append(score, mass);
      if (group.mirrors.length > 1) {
        const layer = document.createElement("select");
        layer.className = "viz-layer-select";
        layer.setAttribute("aria-label", `Library spectrum layer for ${group.name}`);
        group.mirrors.forEach((item, index) => {
          const option = document.createElement("option");
          option.value = String(index);
          option.textContent = `Layer ${index + 1} · score ${fmt(item.score)} · CE ${fmt(item.ce, 1)}`;
          layer.append(option);
        });
        meta.append(layer);
      }
      head.append(title, meta);
      const canvas = document.createElement("canvas");
      canvas._mirror = mirror;
      canvas.setAttribute("aria-label", `Layered mirror spectrum for ${group.name}`);
      const tip = document.createElement("div");
      tip.className = "viz-tooltip is-hidden";
      card.style.position = "relative";
      const select = meta.querySelector("select");
      const showLayer = (index) => {
        canvas._mirror = group.mirrors[index];
        score.textContent = `score ${fmt(canvas._mirror.score)}`;
        mass.textContent = `${fmt(canvas._mirror.libraryMz, 4)} m/z`;
        if (canvas.dataset.live === "true") drawMirror(canvas, canvas._mirror);
      };
      if (select) select.addEventListener("change", () => showLayer(Number(select.value)));
      showLayer(0);
      canvas.addEventListener("pointermove", (event) => mirrorPointer(event, canvas, tip, canvas._mirror));
      canvas.addEventListener("pointerleave", () => hideTip(tip));
      card.append(head, canvas, tip);
      els.mirrorList.append(card);
      mirrorObserver.observe(canvas);
    });
  }

  function drawMirror(canvas, mirror) {
    const { ctx, width, height } = prepareCanvas(canvas, mirrorDpr());
    const p = palette();
    ctx.fillStyle = p.paper;
    ctx.fillRect(0, 0, width, height);
    const box = { left: 54, top: 19, width: Math.max(1, width - 72), height: Math.max(1, height - 48) };
    const baseline = box.top + box.height / 2;
    const maxMz = Math.max(
      mirror.experimentalMz,
      mirror.libraryMz,
      ...mirror.experimental.map((peak) => peak.mz),
      ...mirror.library.map((peak) => peak.mz),
      1,
    ) * 1.03;
    const x = (value) => box.left + (value / maxMz) * box.width;
    const sharedMax = Math.max(
      1,
      ...mirror.experimental.map((peak) => transformedIntensity(peak.intensity)),
      ...mirror.library.map((peak) => transformedIntensity(peak.intensity)),
    );
    ctx.strokeStyle = p.rule;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(box.left, baseline);
    ctx.lineTo(box.left + box.width, baseline);
    ctx.stroke();
    ctx.font = `10px ${p.mono}`;
    ctx.fillStyle = p.ink3;
    ctx.textBaseline = "top";
    ctx.textAlign = "center";
    for (const value of ticks(0, maxMz, 7)) {
      const at = x(value);
      ctx.beginPath();
      ctx.moveTo(at, box.top);
      ctx.lineTo(at, box.top + box.height);
      ctx.stroke();
      ctx.fillText(formatTick(value), at, box.top + box.height + 7);
    }
    const drawSide = (peaks, direction, color) => {
      ctx.strokeStyle = color;
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      for (const peak of peaks) {
        const at = x(peak.mz);
        const end = baseline - direction * (transformedIntensity(peak.intensity) / sharedMax) * (box.height / 2 - 8);
        ctx.moveTo(at, baseline);
        ctx.lineTo(at, end);
      }
      ctx.stroke();
      ctx.fillStyle = color;
      for (const peak of peaks) {
        if (!peak.matched) continue;
        const at = x(peak.mz);
        const end = baseline - direction * (transformedIntensity(peak.intensity) / sharedMax) * (box.height / 2 - 8);
        ctx.beginPath();
        ctx.arc(at, end, 2.6, 0, Math.PI * 2);
        ctx.fill();
      }
    };
    drawSide(mirror.experimental, 1, p.bad);
    drawSide(mirror.library, -1, p.accent);
    ctx.font = `10px ${p.mono}`;
    ctx.textAlign = "left";
    ctx.fillStyle = p.bad;
    ctx.fillText(`experimental · ${fmt(mirror.experimentalMz, 4)} m/z · ${fmt(mirror.rt, 3)} min · CE ${fmt(mirror.ce, 1)}`, box.left, box.top);
    ctx.fillStyle = p.accent;
    ctx.textBaseline = "bottom";
    ctx.fillText(`library · ${fmt(mirror.libraryMz, 4)} m/z`, box.left, box.top + box.height);
    canvas._viz = { box, baseline, x, sharedMax };
  }

  function mirrorPointer(event, canvas, tip, mirror) {
    if (!canvas._viz) return;
    const rect = canvas.getBoundingClientRect();
    const x = event.clientX - rect.left;
    const y = event.clientY - rect.top;
    const list = y < canvas._viz.baseline ? mirror.experimental : mirror.library;
    let best = null;
    let distance = 10;
    for (const peak of list) {
      const next = Math.abs(canvas._viz.x(peak.mz) - x);
      if (next < distance) {
        best = peak;
        distance = next;
      }
    }
    if (!best) {
      hideTip(tip);
      return;
    }
    const side = y < canvas._viz.baseline ? "experimental" : "library";
    showTip(
      tip,
      canvas.parentElement,
      x,
      y + canvas.offsetTop,
      `${side}\n${fmt(best.mz, 4)} m/z · ${formatTick(transformedIntensity(best.intensity))} ${intensityMode === "sqrt" ? "√ intensity" : "intensity"}${best.matched ? " · matched" : ""}`,
    );
  }

  els.sample.addEventListener("change", () => {
    loadToken += 1;
    loadOverview(loadToken);
  });
  els.reload.addEventListener("click", () => loadSession(true));
  els.reset.addEventListener("click", resetMap);
  els.intensityMode.addEventListener("click", (event) => {
    const button = event.target.closest("button");
    if (!button || button.disabled || button.dataset.value === intensityMode) return;
    intensityMode = button.dataset.value === "sqrt" ? "sqrt" : "raw";
    Array.from(els.intensityMode.children).forEach((item) =>
      item.classList.toggle("is-on", item === button)
    );
    try {
      localStorage.setItem("mk-viz-intensity", intensityMode);
    } catch (error) {
      /* storage is optional */
    }
    if (detail) {
      drawSpectrum();
      els.mirrorList.querySelectorAll('canvas[data-live="true"]').forEach((canvas) =>
        drawMirror(canvas, canvas._mirror)
      );
    }
  });

  let detailFrame = 0;
  const requestDetailDraw = () => {
    if (detailFrame) return;
    detailFrame = requestAnimationFrame(() => {
      detailFrame = 0;
      requestMapDraw();
      if (detail) {
        drawXic();
        drawSpectrum();
      }
      els.mirrorList.querySelectorAll('canvas[data-live="true"]').forEach((canvas) => {
        if (canvas._mirror) drawMirror(canvas, canvas._mirror);
      });
    });
  };
  const resize = new ResizeObserver(() => {
    requestDetailDraw();
  });
  [els.map, els.xic, els.spectrum].forEach((canvas) => resize.observe(canvas));

  const redrawForTheme = () => {
    requestMapDraw();
    if (detail) {
      drawXic();
      drawSpectrum();
      els.mirrorList.querySelectorAll('canvas[data-live="true"]').forEach((canvas) => {
        if (canvas._mirror) drawMirror(canvas, canvas._mirror);
      });
    }
  };
  new MutationObserver(redrawForTheme).observe(document.documentElement, {
    attributes: true,
    attributeFilter: ["data-theme"],
  });
  const scheme = window.matchMedia("(prefers-color-scheme: dark)");
  if (scheme.addEventListener) scheme.addEventListener("change", redrawForTheme);
  else if (scheme.addListener) scheme.addListener(redrawForTheme);

  return {
    async setOutput(path, { reload = false } = {}) {
      if (!reload && path === outputDir) return;
      outputDir = path || null;
      session = null;
      overview = null;
      selected = null;
      detail = null;
      clearMirrors();
      els.sample.textContent = "";
      els.sample.disabled = true;
      els.reload.disabled = !outputDir;
      els.reset.disabled = true;
      if (!outputDir) {
        empty("No visualizer data loaded", "Choose a dataset with completed results first.");
        status("Choose a dataset with completed results.");
      } else if (active) {
        await loadSession(true);
      } else {
        empty("Visualizer ready to load", "Open this tab to check and inspect the run data.");
        status("Visualizer data will be checked when this tab opens.");
      }
    },
    async activate() {
      active = true;
      if (outputDir && !session && !analysisRunning) await loadSession();
      else requestMapDraw();
    },
    setRunning(on) {
      analysisRunning = on;
      if (on) {
        loadToken += 1;
        els.sample.disabled = true;
        els.reload.disabled = true;
        els.reset.disabled = true;
        status("Visualizer paused while analysis updates its caches.");
      } else {
        els.sample.disabled = !session || !session.samples.length;
        els.reload.disabled = !outputDir;
        els.reset.disabled = !overview;
      }
    },
  };
}
