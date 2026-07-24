const { invoke } = window.__TAURI__.core;

const load_all = async function () {
  const mzml_Object = await invoke("get_mzml_list");
  const mzml_sel = document.getElementById("mzml_sel");
  mzml_sel.innerHTML = "";
  mzml_sel.options[0] = new Option("--Select--");
  mzml_Object.forEach((x, i) => {
    mzml_sel.options[mzml_sel.options.length] = new Option(x, i);
  });
  const param_t = await invoke("read_param");
  mzml_sel.onchange = function () {
    gen(mzml_Object[this.value], param_t);
  };
  const refreshb = document.getElementById("refreshb");
  refreshb.onclick = function () {
    gen(mzml_Object[mzml_sel.value], param_t);
  };
};

window.onload = async function () {
  load_all();
};

const width = 854;
const height = 480;
const marginTop = 40;
const marginRight = 60;
const marginBottom = 20;
const marginLeft = 40;

const height1 = 180;
const marginTop1 = 20;
const marginRight1 = 30;
const marginBottom1 = 20;
const marginLeft1 = 40;


async function gen(mzml_f, param_t) {
  const { ms1ms2, i_rt, i_mz } = param_t;
  const ms1feat = await invoke("ms1feat", { bn: mzml_f, ms1ms2 });

  const x = d3
    .scaleLinear()
    .domain(d3.extent(ms1feat, (d) => d.rt))
    .range([marginLeft, width - marginRight]);
  const y = d3
    .scaleLinear()
    .domain([ms1feat[0].ms1mz, ms1feat[ms1feat.length - 1].ms1mz])
    .range([height - marginBottom, marginTop]);
  const svg = d3
    .create("svg")
    .style("background", "white")
    .attr("width", width)
    .attr("height", height)
    .style("cursor", "crosshair")
    .style("border", "solid")
    .style("border-radius", "1em");
  svg
    .append("rect")
    .attr("width", width)
    .attr("height", height)
    .style("fill", "none");
  const xAxis = (g, x) => {
    g.call(d3.axisBottom(x));
    g.select(".domain").remove();
    g.attr("font-size", 12);
  };
  const yAxis = (g, y) => {
    g.call(d3.axisLeft(y));
    g.select(".domain").remove();
    g.attr("font-size", 12);
  };
  const gDot = svg
    .append("g")
    .selectAll("circle")
    .data(ms1feat)
    .join("circle")
    .attr("fill", (d) => (d.name_l.length > 0 ? "red" : "none"))
    .attr("stroke", (d) => (d.name_l.length > 0 ? "red" : "black"))
    .attr("r", 3);
  const delaunay = d3.Delaunay.from(
    ms1feat,
    (d) => x(d.rt),
    (d) => y(d.ms1mz),
  );
  svg
    .append("text")
    .attr("x", 4)
    .attr("dominant-baseline", "text-before-edge")
    .text("↑ m/z");
  svg
    .append("text")
    .attr("x", "100%")
    .attr("y", height - marginBottom)
    .attr("text-anchor", "end")
    .attr("dominant-baseline", "text-before-edge")
    .text("RT →");
  const gx = svg
    .append("g")
    .attr("transform", `translate(0,${height - marginBottom})`);
  const gy = svg.append("g").attr("transform", `translate(${marginLeft})`);
  let transform,
    zx,
    zy,
    selpt = [999, 0];
  const zoom = d3
    .zoom()
    .scaleExtent([1, 999])
    .translateExtent([
      [0, 0],
      [width, height],
    ])
    .on("zoom", (e) => {
      tt.style("display", "none");
      transform = e.transform;
      zx = transform.rescaleX(x); //.interpolate(d3.interpolateRound);
      zy = transform.rescaleY(y); //.interpolate(d3.interpolateRound);
      gDot.attr("cx", (d) => zx(d.rt)).attr("cy", (d) => zy(d.ms1mz));
      mark.attr("transform", `translate(${zx(selpt[0])}, ${zy(selpt[1])})`);
      gx.call(xAxis, zx);
      gy.call(yAxis, zy);
    });

  const mark = svg
    .append("circle")
    .attr("opacity", 0.4)
    .attr("r", "1em")
    .attr("fill", "none")
    .attr("stroke-width", "1em")
    .attr("stroke", "black");
  const tt = svg.append("g").style("display", "none");
  tt.append("line")
    .attr("x1", "-100%")
    .attr("x2", "100%")
    .attr("stroke", "black")
    .attr("stroke-opacity", 0.5);
  tt.append("line")
    .attr("y1", "-100%")
    .attr("y2", "100%")
    .attr("stroke", "black")
    .attr("stroke-opacity", 0.5);
  tt.append("rect")
    .attr("width", 200)
    .attr("height", 34)
    .attr("x", 6)
    .attr("y", -40)
    .attr("rx", 4)
    .attr("ry", 4);
  const ttext = tt
    .append("text")
    .attr("y", tt.select("rect").attr("y"))
    .attr("fill", "white");
  const x_pos = 2 + parseFloat(tt.select("rect").attr("x"));
  ttext
    .append("tspan")
    .attr("x", x_pos)
    .attr("dy", 0)
    .attr("dominant-baseline", "text-before-edge");
  ttext
    .append("tspan")
    .attr("x", x_pos)
    .attr("dy", tt.select("rect").attr("height"))
    .attr("dominant-baseline", "text-after-edge");


  const svg_spec = new Array(40)
    .fill(null)
    .map(() => d3.create("svg").style("background", "white").attr("height", 0));

  svg
    .call(zoom)
    .call(zoom.transform, d3.zoomIdentity)
    .on("pointermove", (event) => {
      const p = transform.invert(d3.pointer(event));
      const i = delaunay.find(...p);
      const ms1feat_i = ms1feat[i];
      tt.style("display", null).attr(
        "transform",
        `translate(${zx(ms1feat_i.rt)}, ${zy(ms1feat_i.ms1mz)})`,
      );
      ttext
        .select("tspan")
        .text(
          `${d3.format(".4f")(ms1feat_i.ms1mz)}, ${d3.format(".2f")(ms1feat_i.rt)}`,
        );
      ttext
        .select("tspan:nth-of-type(2)")
        .text(ms1feat_i.name_l.join(", ").slice(0, 30));
    })
    .on("pointerleave", () => tt.style("display", "none"))
    .on("click", async () => {
      const p = transform.invert(d3.pointer(event));
      const i = delaunay.find(...p);
      const ms1feat_i = ms1feat[i];
      selpt = [ms1feat_i.rt, ms1feat_i.ms1mz];
      mark
        .transition()
        .attr("transform", `translate(${zx(selpt[0])}, ${zy(selpt[1])})`);
      const rtwid = 0.5;
      const msms = await invoke("get_spec", {
        bn: mzml_f,
        ms1mz: ms1feat_i.ms1mz,
        ms1rt: ms1feat_i.rt,
        rtwid,
        ms1ms2,
      });
      const chrom = await invoke("get_ms1", {
        bn: mzml_f,
        ms1mz: ms1feat_i.ms1mz,
        ms1rt: ms1feat_i.rt,
        rtwid,
        iMz: i_mz,
      });
      for (const spec of svg_spec) {
        spec.selectAll("svg > *").remove();
        spec.attr("width", width).attr("height", 0).style("border", "none");
      }
      print_xic(svg_spec, msms, chrom, ms1feat_i, i_rt);
      const mirror = await invoke("get_mirror", {
        bn: mzml_f,
        ms1mz: ms1feat_i.ms1mz,
        ms1rt: ms1feat_i.rt,
      });
      mirror.forEach((msms_, i) => {
        print_mirror(svg_spec[2 + i], msms_);
      });
    });

  const container = document.getElementById("container");
  container.innerHTML = "";
  container.append(svg.node());
  for (const spec of svg_spec) {
    container.append(spec.node());
  }
}
function print_xic(svg_spec, msms, chrom, ms1feat_i, i_rt) {
  const svg = svg_spec[0];
  svg
    .attr("height", height1)
    .style("border", "solid")
    .style("border-radius", "1em");

  svg.selectAll("svg > *").remove();
  svg
    .append("text")
    .attr("x", "50%")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge")
    .attr("font-weight", "bold")
    .text(
      `XIC @ ${d3.format(".4f")(ms1feat_i.ms1mz)}m/z, shape: ${d3.format(".2f")(ms1feat_i.shape)}, SN: ${d3.format(".2f")(ms1feat_i.smooth)}`,
    );
  const x = d3
    .scaleLinear()
    .domain([chrom[0][0], chrom[chrom.length - 1][0]])
    .range([marginLeft1, width - marginRight1]);
  const y = d3
    .scaleLinear()
    .domain([0, 1.05 * d3.max(chrom, (d) => d[1])])
    .range([height1 - marginBottom1, marginTop1]);
  svg
    .append("g")
    .attr("transform", `translate(0,${height1 - marginBottom1})`)
    .call(d3.axisBottom(x))
    .call((g) => g.select(".domain").attr("opacity", 0.5))
    .call((g) => g.attr("font-size", 12));
  svg
    .append("g")
    .attr("transform", `translate(${marginLeft1})`)
    .call(d3.axisLeft(y).ticks(2, "s"))
    .call((g) => g.select(".domain").remove())
    .call((g) => g.attr("font-size", 12));
  svg
    .append("g")
    .selectAll("circle")
    .data(chrom)
    .join("circle")
    .attr("cx", (d) => x(d[0]))
    .attr("cy", (d) => y(d[1]))
    .attr("r", 2);
  const gLine = svg
    .append("g")
    .attr("opacity", 0.5)
    .attr("stroke", "red")
    .selectAll("line")
    .data(msms)
    .join("line")
    .attr("x1", (d) => x(d.rt))
    .attr("x2", (d) => x(d.rt))
    .attr("y1", y.range()[0])
    .attr("y2", y.range()[1]);
  const tt = svg.append("g").style("display", "none");
  tt.append("rect")
    .attr("width", 90)
    .attr("height", marginBottom1)
    .attr("x", -tt.select("rect").attr("width") / 2)
    .attr("y", y.range()[0])
    .attr("rx", 4)
    .attr("ry", 4);
  tt.append("text")
    .attr("fill", "white")
    .attr("y", y.range()[0])
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge");

  const int_wid = ms1feat_i.sc * i_rt;
  const int_begin = x(ms1feat_i.rt - int_wid);
  svg
    .append("rect")
    .attr("width", x(ms1feat_i.rt + int_wid) - int_begin)
    .attr("height", y.range()[0] - y.range()[1])
    .attr("x", int_begin)
    .attr("y", y.range()[1])
    .attr("opacity", 0.1);

  const line = d3
    .line()
    .x((d) => x(d[0]))
    .y((d) => y(d[1]));
  svg
    .append("path")
    .attr("d", line(chrom))
    .attr("fill", "none")
    .attr("stroke", "black");
  const bisect = d3.bisector((d) => d.rt).center;
  svg_spec[1]
    .attr("width", width)
    .attr("height", height1)
    .style("border", "solid")
    .style("border-radius", "1em");
  svg
    .on("pointermove", (event) => {
      const i = bisect(msms, x.invert(d3.pointer(event)[0]));
      gLine.attr("stroke-width", (d, ii) => (ii == i ? 3 : 1));
      tt.style("display", null).attr(
        "transform",
        `translate(${x(msms[i].rt)})`,
      );
      tt.select("text").text(
        d3.format(".4f")(msms[i].ms1mz - ms1feat_i.ms1mz) + "m/z",
      );
      print_spec(svg_spec[1], msms[i]);
    })
    .on("pointerleave", () => tt.style("display", "none"));
}
function print_mirror(svg, mirror) {
  svg
    .attr("height", height1)
    .style("border", "solid")
    .style("border-radius", "1em");
  svg
    .append("text")
    .attr("x", "50%")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge")
    .attr("font-weight", "bold")
    .attr("fill", "red")
    .text(
      `MS/MS @ ${d3.format(".4f")(mirror.specmz)}m/z, ${d3.format(".3f")(mirror.specrt)}min, CE: ${d3.format(".1f")(mirror.ce)}`,
    );
  svg
    .append("text")
    .attr("x", "99%")
    .attr("text-anchor", "end")
    .attr("dominant-baseline", "text-before-edge")
    .attr("font-weight", "bold")
    .text(`score: ${d3.format(".2f")(mirror.dotp)}`);
  svg
    .append("text")
    .attr("x", "50%")
    .attr("y", "100%")
    .attr("text-anchor", "middle")
    .attr("font-weight", "bold")
    .attr("dominant-baseline", "ideographic")
    .attr("fill", "blue")
    .text(mirror.name + ` (${d3.format(".4f")(mirror.lib_mass)}m/z)`);
  const x = d3
    .scaleLinear()
    .domain([0, mirror.ms1mz + 0.3])
    .range([marginLeft1, width - marginRight1]);
  const y = d3
    .scaleLinear()
    .domain([-1, 0, mirror.exp_mz_i_l[0][1]])
    .range([height1 - marginBottom1 - 5, height1 / 2, marginTop1 + 5]);
  const all_p = mirror.exp_mz_i_l.concat(
    mirror.lib_mz_i_l.map((d) => [d[0], -d[1]]),
  );
  svg
    .append("g")
    .attr("stroke-width", 2)
    .selectAll("line")
    .data(all_p)
    .join("line")
    .attr("x1", (d) => x(d[0]))
    .attr("y1", y(0))
    .attr("x2", (d) => x(d[0]))
    .attr("y2", (d) => y(d[1]))
    .attr("stroke", (d) => (d[1] > 0 ? "red" : "blue"));
  const ma_p = mirror.m_exp_mz_i_l.concat(
    mirror.m_lib_mz_i_l.map((d) => [d[0], -d[1]]),
  );
  svg
    .append("g")
    .selectAll("circle")
    .data(ma_p)
    .join("circle")
    .attr("cx", (d) => x(d[0]))
    .attr("cy", (d) => y(d[1]))
    .attr("fill", (d) => (d[1] > 0 ? "red" : "blue"))
    .attr("r", 3);
  const bisect = d3.bisector((d) => d[0]).center;
  const tt = svg.append("g").style("display", "none");
  const tt1 = tt.append("g");
  tt1
    .append("rect")
    .attr("width", 70)
    .attr("height", 20)
    .attr("fill", "red")
    .attr("x", -tt1.select("rect").attr("width") / 2)
    .attr("y", -24)
    .attr("rx", 4)
    .attr("ry", 4);
  tt1
    .append("text")
    .attr("y", -14)
    .attr("fill", "white")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "central");
  const tt2 = tt.append("g");
  tt2
    .append("rect")
    .attr("width", 70)
    .attr("height", 20)
    .attr("fill", "blue")
    .attr("x", -tt2.select("rect").attr("width") / 2)
    .attr("y", 4)
    .attr("rx", 4)
    .attr("ry", 4);
  tt2
    .append("text")
    .attr("y", 14)
    .attr("fill", "white")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "central");
  svg
    .on("pointermove", (event) => {
      const i = bisect(mirror.m_exp_mz_i_l, x.invert(d3.pointer(event)[0]));
      tt.style("display", null);
      tt1.attr(
        "transform",
        `translate(${x(mirror.m_exp_mz_i_l[i][0])}, ${y(mirror.m_exp_mz_i_l[i][1])})`,
      );
      tt1.select("text").text(d3.format(".3f")(mirror.m_exp_mz_i_l[i][0]));
      tt2.attr(
        "transform",
        `translate(${x(mirror.m_lib_mz_i_l[i][0])}, ${y(-mirror.m_lib_mz_i_l[i][1])})`,
      );
      tt2.select("text").text(d3.format(".3f")(mirror.m_lib_mz_i_l[i][0]));
    })
    .on("pointerleave", () => tt.style("display", "none"));
  svg
    .append("line")
    .attr("y1", y(0))
    .attr("x2", "100%")
    .attr("y2", y(0))
    .attr("opacity", 0.5)
    .attr("stroke", "black");
}
function print_spec(svg, { mz_i_l, ms1mz, rt, ce }) {
  svg
    .attr("width", width)
    .attr("height", height1)
    .style("border", "solid")
    .style("border-radius", "1em");
  svg.selectAll("svg > *").remove();
  svg
    .append("rect")
    .attr("width", width)
    .attr("height", height1)
    .style("fill", "none");
  svg
    .append("text")
    .attr("x", "50%")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "text-before-edge")
    .attr("font-weight", "bold")
    .text(
      `MS/MS @ ${d3.format(".4f")(ms1mz)}m/z, ${d3.format(".3f")(rt)}min, CE: ${d3.format(".1f")(ce)}`,
    );
  const x = d3
    .scaleLinear()
    .domain([0, mz_i_l[mz_i_l.length - 1][0]])
    .range([marginLeft1, width - marginRight1]);
  const y = d3.scaleLinear().range([height1 - marginBottom1, marginTop1]);
  const gLine = svg
    .append("g")
    .attr("stroke-width", 2)
    .attr("stroke", "black")
    .selectAll("line")
    .data(mz_i_l)
    .join("line")
    .attr("y1", y.range()[0]);
  const mark = svg
    .append("polygon")
    .attr("points", "0,0 -11.547,20 11.547,20")
    .attr("opacity", 0.5)
    .attr("fill", "red");
  const xAxis = (g, x) => {
    g.call(d3.axisBottom(x));
    g.select(".domain").attr("opacity", 0.5);
    g.attr("font-size", 12);
  };
  const yAxis = (g, y) => {
    g.call(d3.axisLeft(y).ticks(2, "s"));
    g.select(".domain").remove();
    g.attr("font-size", 12);
  };
  const gx = svg
    .append("g")
    .attr("transform", `translate(0,${height1 - marginBottom1})`);
  const gy = svg.append("g").attr("transform", `translate(${marginLeft1})`);
  let xz;
  const bi = d3.bisector((d) => d[0]);
  const zoom = d3
    .zoom()
    .translateExtent([
      [0, -Infinity],
      [width, Infinity],
    ])
    .scaleExtent([1, 999])
    .on("zoom", (e) => {
      tt.style("display", "none");
      xz = e.transform.rescaleX(x);
      const pos0 = bi.right(mz_i_l, xz.invert(marginLeft1));
      const pos1 = bi.right(mz_i_l, xz.invert(width - marginRight1));
      const max_y = d3.max(mz_i_l.slice(pos0, pos1), (d) => d[1]);
      y.domain([0, max_y === undefined ? 1 : 1.05 * max_y]);
      mark.attr("transform", `translate(${xz(ms1mz)}, ${y.range()[0]})`);
      gLine
        .attr("x1", (d) => xz(d[0]))
        .attr("x2", (d) => xz(d[0]))
        .attr("y2", (d) => y(d[1]));
      gx.call(xAxis, xz);
      gy.call(yAxis, y);
      delaunay = d3.Delaunay.from(
        mz_i_l,
        (d) => xz(d[0]),
        (d) => y(d[1]),
      );
    });
  const tt = svg.append("g").style("display", "none");
  tt.append("rect")
    .attr("width", 70)
    .attr("height", 20)
    .attr("x", -tt.select("rect").attr("width") / 2)
    .attr("y", -24)
    .attr("rx", 4)
    .attr("ry", 4);
  tt.append("text")
    .attr("y", -14)
    .attr("fill", "white")
    .attr("text-anchor", "middle")
    .attr("dominant-baseline", "central");
  tt.append("circle").attr("r", 3);
  let delaunay;
  svg
    .call(zoom)
    .call(zoom.transform, d3.zoomIdentity)
    .on("pointermove", (event) => {
      const i = delaunay.find(...d3.pointer(event));
      tt.style("display", null).attr(
        "transform",
        `translate(${xz(mz_i_l[i][0])}, ${y(mz_i_l[i][1])})`,
      );
      tt.select("text").text(d3.format(".3f")(mz_i_l[i][0]));
    })
    .on("pointerleave", () => tt.style("display", "none"));
}
