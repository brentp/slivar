const ROC_VAR = "GQ"
var sliders = { samples: {} }
var plots = { hists: {}, rocs: {} }
var stats = {}
var sorted = {}
const colors = ["#4e79a7", "#f28e2c", "#e15759", "#76b7b2", "#59a14f", "#edc949", "#af7aa1", "#ff9da7", "#9c755f", "#bab0ab"]
const primary_color = "#2F3C48"
const main_roc_layout = {
    height: 400,
    margin: { t: 10, r: 0, b: 0, l: 60 },
    xaxis: { automargin: true, title: { text: "Mendelian Violations", standoff: 100 }, rangemode: 'tozero', },
    yaxis: { automargin: true, title: { text: "Transmitted Variants", standoff: 10 }, rangemode: 'tozero', },
    font: { size: 15 }
}
const aux_plot_height = 250
const histogram_layout = {
    barmode: "overlay",
    xaxis: { automargin: true, },
    // xaxis: { automargin: true, rangeslider: {} },
    yaxis: { title: { text: "Count", standoff: 20 }, autorange: true, automargin: true, type: "log", fixedrange: true },
    height: aux_plot_height,
    margin: { t: 10, b: 20, r: 0, l: 50, pad: 0 },
    legend: {
        xanchor: "right",
        yanchor: "top",
        y: 1.1,
        x: 1,
        orientation: "h",
    },
    // dragmode: 'select'
}
const roc_layout = {
    xaxis: { automargin: true, title: { text: "Violations", standoff: 20 } },
    yaxis: { automargin: true, title: { text: "Transmitted", standoff: 20 } },
    height: aux_plot_height,
    margin: { t: 10, b: 20, r: 10, l: 30, pad: 0 },
    legend: {
        xanchor: "right",
        yanchor: "top",
        y: 1.2,
        x: 1,
        orientation: "h",
    },
}
const box_layout = {
    barmode: 'group',
    height: aux_plot_height,
    xaxis: { automargin: true, },
    yaxis: { title: { standoff: 20 }, autorange: true, automargin: true, fixedrange: true },
    margin: { t: 10, b: 20, r: 0, l: 80, pad: 0 },
    legend: {
        xanchor: "right",
        yanchor: "top",
        y: 1.1,
        x: 1,
        orientation: "h",
    },
}


function set_stats() {
    stats.n_violations = 0
    for (let i = 0; i < variant_infos.violations.length; i++) {
        stats.n_violations += variant_infos.violations[i];
    }
    stats.n_transmitted = variant_infos.violations.length - stats.n_violations;
}

function update_stats(idxs) {
    let fpr = (stats.n_violations / (stats.n_violations + stats.n_transmitted)).toFixed(3)
    if (stats.filtered == undefined) {
        stats.filtered = {}
    }
    stats.filtered.n_violations = 0
    stats.filtered.n_variants = idxs.size
    stats.filtered.n_transmitted = 0
    idxs.forEach(function (e) {
        stats.filtered.n_violations += variant_infos.violations[e]
    })
    stats.filtered.n_transmitted = idxs.size - stats.filtered.n_violations
    // no filters have been applied
    if (isNaN(stats.filtered.n_transmitted)) {
        jQuery('#metrics-transmitted').html(`${stats.n_transmitted}`)
        jQuery('#metrics-violations').html(`${stats.n_violations}`)
        jQuery('#metrics-fpr').html(fpr)
    } else {
        // will need to know when all filters have been removed
        let f_fpr
        if (stats.filtered.n_violations + stats.filtered.n_transmitted > 0) {
            f_fpr = (stats.filtered.n_violations / (stats.filtered.n_violations + stats.filtered.n_transmitted)).toFixed(3)
        } else {
            f_fpr = 0
        }
        jQuery('#metrics-transmitted').html(`${stats.filtered.n_transmitted} <span class="text-info">out of ${stats.n_transmitted} (${(100 * stats.filtered.n_transmitted / stats.n_transmitted).toFixed(2)}%)</span>`)
        jQuery('#metrics-violations').html(`${stats.filtered.n_violations} <span class="text-info">out of ${stats.n_violations} (${(100 * stats.filtered.n_violations / stats.n_violations).toFixed(2)}%)</span>`)
        jQuery('#metrics-fpr').html(`${f_fpr} <span class="text-info">from ${fpr}</span>`)
    }
}

function iNumb(opts) {
    return {
        to: function (a) {
            return parseInt(Math.round(a)).toString()
        },
        from: function (txt) {
            return parseInt(txt)
        },
    }
}

function wNumb(opts) {
    return {
        to: function (a) {
            return a.toFixed(opts.decimals).replace(/([0-9]+(\.[0-9]+[1-9])?)(\.?0+$)/, '$1')
        },
        from: function (txt) {
            return parseFloat(txt);
        }
    };
}

function roc(values, violations, invert, name, filters_to_keep, idxs) {
    let filters = variant_infos.filters

    //console.time("deco-sort " + name)
    let key = name + "__" + invert
    if (sorted[key] == undefined) {
        sorted[key] = values.map(function (val, i) {
            if (invert) { val = -val };
            return [val, violations[i]]
        })
        sorted[key].sort(function (a, b) {
            return a[0] - b[0]
        })
    }
    var Aorig = sorted[key];
    var Afilt = Aorig.filter(function (val, i) {
        return (idxs.size == 0 || idxs.has(i)) && filters_to_keep.includes(filters[i])
    })
    //console.timeEnd("deco-sort " + name)
    let N = Aorig.length
    var txt = "INFO." + name + (invert ? ">" : "<") + " "
    // we can't just iterate over A, we also have to track the change in
    // cutoff so we don't draw a smooth curve when the cutoff value hasn't
    // changed.
    result = [
        { x: [0], y: [0], text: [""], name: "Original", marker: { color: colors[2] }, line: { dash: 'solid', width: 2 } },
        { x: [0], y: [0], text: [""], name: "Filtered", marker: { color: colors[3] }, line: { dash: 'solid', width: 2 } }
    ]
    if (N == 0) {
        return result
    }
    //console.time("prepare traces " + name) ;

    [Aorig, Afilt].forEach(function (A, Ai) {
        var last_val = A[0][0] - 1;
        var tps = 0
        var fps = 0
        A.forEach(function (pair, i) {
            //if(pair[1]) { fps += 1; }
            //else {tps += 1; }
            fps += pair[1] // pair[1] indicates that it's a violation
            tps += (1 - pair[1])
            if (pair[0] != last_val) {
                let val = pair[0]
                last_val = val;

                // TODO: scale by A.length
                // but might be good to leave as numbers for now.
                result[Ai].x.push(fps)
                result[Ai].y.push(tps)
                result[Ai].text.push(`${txt} ${val}. FDR: ${(fps / (tps + fps)).toFixed(3)}`)
            }
        })
        result[Ai].x.push(fps)
        result[Ai].y.push(tps)
        result[Ai].text.push(`${txt} ${A[A.length - 1][0]} FDR: ${(fps / (tps + fps)).toFixed(3)}`)
    })
    //console.timeEnd("prepare traces " + name)

    return result
}

function arr_min_max(arr) {
    var min = arr[0];
    var max = arr[0]
    for (i = 1; i < arr.length; i++) {
        if (arr[i] < min) { min = arr[i]; }
        if (arr[i] > max) { max = arr[i]; }
    }
    if (Math.abs(parseInt(min) - min) > 0.001) {
        min -= 0.004;
        max += 0.004;
    }
    return [min, max]
}

function add_bool(values, name, idxs) {
    jQuery('#variant-filter-nav').append(`
        <a class="nav-link ml-3 my-1" href="#${name}-wrapper">${name}</a>
    `)
    jQuery('#variant-sliders').append(`
        <div class="container-fluid pb-5 pt-3">
            <a id=${name}-wrapper></a>
            <label for="${name}-bar-plot"><h3>${name}</h3></label>
            <!-- Plot area -->
            <div class="row border-top">
                <div class="col-12 p-1" id=${name}-bar-plot></div>
            </div>
            <!-- Boolean selection -->
            <div class="row pb-3 border-bottom">
                <div class="col-1 pr-0">
                    <span class="selection-label">Include:</span>
                </div>
                <div class="col-11">
                    <div class="btn-group-toggle" data-toggle="buttons">
                        <label class="btn btn-outline-primary m-1 p-1 active">
                            <input type="checkbox" autocomplete="off" class="slivar-checkbox slivar-bool-checkbox slivar-changer" checked name=${name}-yes id=${name}-yes>
                                Yes (True)
                            </input>
                        </label>
                        <label class="btn btn-outline-primary m-1 p-1 active">
                            <input type="checkbox" autocomplete="off" class="slivar-checkbox slivar-bool-checkbox slivar-changer" checked name=${name}-no id=${name}-no>
                                No (False)
                            </input>
                        </label>
                    </div>
                </div>
            </div>
        </div>
    `)
    plot_bool(values, name, variant_infos.violations, idxs)
}

function add_slider(values, name, label, is_fmt_field, idxs, violations) {
    var prefix = is_fmt_field ? "sample-" : ""

    // add the filter div
    if (is_fmt_field) {
        jQuery('#sample-filter-nav').append(`
            <a class="nav-link ml-3 my-1" href="#${prefix}${name}-wrapper">${label}</a>
        `)
        jQuery('#sample-sliders').append(`
            <div class="container-fluid pb-5 pt-3">
                <a id=${prefix}${name}-wrapper></a>
                <label for="${prefix}${name}-slider"><h3>${label}</h3></label>
                <!-- Plot area -->
                <div class="row border-top">
                    <div class="col-12 p-1" id=${prefix}${name}-hist-plot></div>
                </div>
                <!-- Range selection -->
                <div class="row pb-5 border-bottom">
                    <div class="col-1 pr-0">
                        <span class="selection-label">Include:</span>
                    </div>
                    <div class="col-11">
                        <div id=${prefix}${name}-slider name=${prefix}${name}-slider></div>
                    </div>
                </div>
            </div>
        `)
    } else {
        jQuery('#variant-filter-nav').append(`
            <a class="nav-link ml-3 my-1" href="#${prefix}${name}-wrapper">${label}</a>
        `)
        jQuery('#variant-sliders').append(`
            <div class="container-fluid pb-5 pt-3">
                <a id=${prefix}${name}-wrapper></a>
                <label for="${prefix}${name}-slider"><h3>${label}</h3></label>
                <!-- Plot area -->
                <div class="row border-top">
                    <div class="col-8 p-1" id=${prefix}${name}-hist-plot></div>
                    <div class="col-4 p-1" id=${prefix}${name}-roc-plot></div>
                </div>
                <!-- Range selection -->
                <div class="row pb-5 border-bottom">
                    <div class="col-1 pr-0">
                        <span class="selection-label">Include:</span>
                    </div>
                    <div class="col-7">
                        <div id=${prefix}${name}-slider name=${prefix}${name}-slider></div>
                    </div>
                    <div class="col-4 text-right">
                        <div class="btn-group-toggle" data-toggle="buttons">
                            <label class="btn btn-sm btn-outline-primary">
                                <input type="checkbox" autocomplete="off" class="slivar-inverter" id="${prefix}${name}-flip" name="${prefix}${name}">
                                    Flip Values
                                </input>
                            </label>
                        </div>
                    </div>
                </div>

            </div>
        `)
    }

    var [vlmin, vlmax] = arr_min_max(values); // can't use apply or destructring as we run out of stack
    var fmtF = name == "variant-length" ? iNumb : wNumb

    sliders[`${prefix}${name}`] = noUiSlider.create(document.getElementById(`${prefix}${name}-slider`), {
        tooltips: [fmtF({ decimals: 3 }), fmtF({ decimals: 3 })],
        start: [vlmin, vlmax],
        connect: true,
        range: { min: vlmin, max: vlmax },
        pips: {
            mode: 'range',
            density: 5,
        }
    })
    // add min and max values to the slider object so we can tell if it has
    // been modified
    sliders[`${prefix}${name}`].min = vlmin
    sliders[`${prefix}${name}`].max = vlmax
    sliders[`${prefix}${name}`].on('change', function () {
        let idxs = new Set(get_passing_idxs())
        main_plot(idxs)
        update_info_plots(idxs)
    })

    plot_field(values, name, label, is_fmt_field, violations, idxs)
}

function get_sample_bounds() {
    var sample_bounds = {}
    for (k in sliders) {
        if (k == "samples") { continue; }
        if (!k.startsWith("sample-")) {
            continue
        }
        let tmp = sliders[k].get()
        if (Math.abs(parseFloat(tmp[0]) - sliders[k].min) < 0.005 && Math.abs(parseFloat(tmp[1]) - sliders[k].max) < 0.005) {
            continue;
        }
        // remove sample- prefix
        sample_bounds[k.substring(7)] = [parseFloat(tmp[0]), parseFloat(tmp[1])]
    }
    return sample_bounds
}

function get_excluded_sample_idxs(sample_bounds) {
    // if a variant is removed in any trio, then we remove it from all
    // consideration
    var result = []
    // console.log(sample_bounds)
    for (field in sample_bounds) {
        let [lo, hi] = sample_bounds[field]
        for (var k = 0; k < trios.length; k++) {
            var trio = trios[k];
            // console.log(trio)
            var values = trio.tbl[field]
            for (var i = 0; i < values.length; i++) {
                if ((values[i] < lo) || (values[i] > hi)) {
                    result.push(i)
                }
            }
        }
    }
    return new Set(result)
}

function get_passing_idxs() {
    var info_bounds = {}
    var info_bools = {}
    let sample_bounds = get_sample_bounds()

    for (k in sliders) {
        if (k == "samples") { continue; }
        if (k.startsWith("sample-")) {
            continue
        } else {
            var tmp = sliders[k].get()
            if (Math.abs(parseFloat(tmp[0]) - sliders[k].min) < 0.005 && Math.abs(parseFloat(tmp[1]) - sliders[k].max) < 0.005) {
                continue
            }
            info_bounds[k] = [parseFloat(tmp[0]), parseFloat(tmp[1])]
        }
    }

    jQuery('.slivar-bool-checkbox').each(function () {
        var n = this.name
        var i = -1;
        if (n.endsWith("-yes")) {
            n = n.substring(0, n.length - 4)
            i = 1
        } else if (n.endsWith("-no")) {
            n = n.substring(0, n.length - 3)
            i = 0
        } else {
            alert("unknown bool type")
        }
        if (!(n in info_bools)) { info_bools[n] = [false, false]; }

        // so now, if info_bools[$x][0] is true, then we include things without
        // the flag. if [1] is true, we include things with the flag.
        info_bools[n][i] = this.checked

    })
    // console.log(info_bools)


    var filters_to_keep = jQuery('.slivar-filters:checked').map(function (v) { return this.name }).toArray();

    var variant_lengths_ranges = info_bounds['variant-length']
    delete info_bounds['variant-length'];

    var excluded_sample_idxs = get_excluded_sample_idxs(sample_bounds)
    let filters = variant_infos.filters;
    var idxs = []
    for (var i = 0; i < variant_infos.variant_lengths.length; i++) {
        var vl = variant_infos.variant_lengths[i];
        if (variant_lengths_ranges != undefined && (vl < variant_lengths_ranges[0] || vl > variant_lengths_ranges[1])) { continue; }
        if (!filters_to_keep.includes(filters[i])) { continue; }

        var skip = false;
        for (k in info_bounds) {
            var vf = variant_infos.float_tbl[k][i];
            if (vf < info_bounds[k][0] || vf > info_bounds[k][1]) {
                skip = true;
                break;
            }
        }
        if (skip) { continue }
        for (k in info_bools) {
            var fl = variant_infos.bool_tbl[k][i];
            if (!info_bools[k][Number(fl)]) {
                skip = true;
                break;
            }
        }
        if (excluded_sample_idxs.has(i)) {
            continue
        }
        if (!skip) { idxs.push(i) }
    }
    return idxs
}

function adjust_with_order(arr, order) {
    var result = new Array(arr.length)

    for (var i = 0; i < arr.length; i++) {
        result[i] = arr[order[i]];
    }
    return result;
}

function sort_trio(trio) {
    // on page-load, we sort each trio so that the kid's lowest ABs
    // are first so we can draw the ROC curve
    var order = new Array(trio.violations.length)
    for (var i = 0; i < trio.violations.length; i++) { order[i] = i; }

    var roc_var = trio.tbl[ROC_VAR];
    //order.sort(function(a, b) { return roc_var[a] - roc_var[b]; });
    order.sort(function (a, b) { return (roc_var[b] - roc_var[a]); });

    for (k in trio.tbl) {
        trio.tbl[k] = adjust_with_order(trio.tbl[k], order)
    }
    trio.violations = adjust_with_order(trio.violations, order);
    trio.variant_idxs = adjust_with_order(trio.variant_idxs, order);
}


function sample_skippable(trio, sample_bounds, i) {
    for (k in sample_bounds) {
        let val = trio.tbl[k][i]
        if (val < sample_bounds[k][0] || val > sample_bounds[k][1]) { return true; }
    }
    return false;
}

function main_plot(idxs) {
    // main-roc-plot
    console.time('main-plot')
    var traces = [];
    // we only plot points that pass the INFO filters
    if (idxs == undefined) {
        idxs = new Set(get_passing_idxs())
    }
    let sample_bounds = get_sample_bounds()
    trios.forEach(function (trio, index) {
        var tps = 0
        var fps = 0
        var last_val = trio.tbl[ROC_VAR] - 1;
        var trace = { x: [], y: [], text: [], name: trio.sample_id } //, line: { color: colors[index % colors.length]} }

        for (let i = 0; i < trio.variant_idxs.length; i++) {
            if (!idxs.has(trio.variant_idxs[i])) { continue; }

            if (sample_skippable(trio, sample_bounds, i)) { continue; }
            var vio = trio.violations[i] //== 0 && trio.mom_alts[i] == 0 && trio.kid_alts[i] == 1;
            if (vio) { fps += 1; }
            else { tps += 1; }
            if (trio.tbl[ROC_VAR][i] == last_val) { continue; }
            last_val = trio.tbl[ROC_VAR][i];

            trace.x.push(fps);
            trace.y.push(tps);
            trace.text.push(`${ROC_VAR} > ${last_val} fpr: ${(fps / (fps + tps)).toFixed(3)}`)

        }
        trace.x.push(fps);
        trace.y.push(tps);
        trace.text.push(`${ROC_VAR} > ${last_val} fpr: ${(fps / (fps + tps)).toFixed(3)}`)
        traces.push(trace)
    })
    console.timeEnd('main-plot')

    Plotly.newPlot('main-roc-plot', traces, main_roc_layout)
}

function plot_bool(values, name, vios, idxs) {
    var filters_to_keep = jQuery('.slivar-filters:checked').map(function (v) { return this.name }).toArray();
    // convert bool to index (false==0) and increment counter for passing
    // variants
    var inh_counts = [0, 0]
    var vio_counts = [0, 0]

    var filt_inh_counts = [0, 0]
    var filt_vio_counts = [0, 0]
    let filters = variant_infos.filters;
    for (var k = 0; k < values.length; k++) {
        if (vios[k]) {
            vio_counts[Number(values[k])]++
        } else {
            inh_counts[Number(values[k])]++
        }
        if (idxs.size > 0 && !idxs.has(k)) { continue; }
        if (!filters_to_keep.includes(filters[k])) { continue; }
        if (vios[k]) {
            filt_vio_counts[Number(values[k])]++
        } else {
            filt_inh_counts[Number(values[k])]++
        }
    }
    var Y = jQuery(`#${name}-yes`).is(":checked")
    var N = jQuery(`#${name}-no`).is(":checked")

    // unfilt, filt
    var vio = [0, 0]
    var tra = [0, 0]
    if (Y) {
        vio[0] += vio_counts[1]
        tra[0] += inh_counts[1]

        vio[1] += filt_vio_counts[1]
        tra[1] += filt_inh_counts[1]
    }
    if (N) {
        vio[0] += vio_counts[0]
        tra[0] += inh_counts[0]

        vio[1] += filt_vio_counts[0]
        tra[1] += filt_inh_counts[0]
    }

    let xlabels = ["Unfiltered", "Filtered"]
    var traces = [
        { x: tra, y: xlabels, name: "Transmitted", type: "bar", orientation: 'h', marker: { color: colors[0] }, },
        { x: vio, y: xlabels, name: "Violations", type: "bar", orientation: 'h', marker: { color: colors[1] }, },
    ]
    Plotly.newPlot(`${name}-bar-plot`, traces, box_layout)
}

function plot_field(values, name, label, is_fmt_field, vios, idxs) {
    // update_only is used only when is_fmt_field to add samples to an existing
    // fmt roc plot
    var [vlmin, vlmax] = arr_min_max(values);
    var filters_to_keep = jQuery('.slivar-filters:checked').map(function (v) { return this.name }).toArray();

    var inh_vals = []
    var vio_vals = []
    let filters = variant_infos.filters;
    //console.time("plot_field filtering:" + name)
    for (var k = 0; k < values.length; k++) {
        if (idxs.size > 0 && !idxs.has(k)) { continue; }
        if (!filters_to_keep.includes(filters[k])) { continue; }
        if (vios[k]) {
            vio_vals.push(values[k])
        } else {
            inh_vals.push(values[k])
        }
    }

    //console.timeEnd("plot_field filtering:" + name)
    var prefix = is_fmt_field ? "sample-" : "";
    var do_flip = jQuery(`#${prefix}${name}-flip`).is(":checked")
    var bin_size = (vlmax - vlmin) / 100;
    var roc_trs;
    if (!is_fmt_field) {
        roc_trs = roc(values, vios, do_flip, name, filters_to_keep, idxs);
    }

    var traces = [
        { type: 'histogram', x: inh_vals, name: "Transmitted", autobinx: false, xbins: { size: bin_size }, histnorm: "count", marker: { color: colors[0] } }, //, opacity: 0.5 },
        { type: 'histogram', x: vio_vals, name: "Violations", autobinx: false, xbins: { size: bin_size }, histnorm: "count", marker: { color: colors[1] } }, //, opacity: 0.5 },
    ];

    plots.hists[`${prefix}${name}`] = Plotly.newPlot(`${prefix}${name}-hist-plot`, traces, histogram_layout, { displayModeBar: false });
    if (!is_fmt_field) {
        plots.rocs[`${prefix}${name}`] = Plotly.newPlot(`${prefix}${name}-roc-plot`, roc_trs, roc_layout, { displayModeBar: false });
    }
}

function update_info_plots(idxs) {
    if (idxs == undefined) {
        idxs = new Set(get_passing_idxs())
    } else {
        console.log("got idxs", idxs)
    }

    let prefix = "";
    update_stats(idxs)
    plot_field(variant_infos.variant_lengths, "variant-length", "variant-lengths", false, variant_infos.violations, idxs);
    for (k in variant_infos.float_tbl) {
        plot_field(variant_infos.float_tbl[k], k, k, false, variant_infos.violations, idxs);
    }
    for (k in variant_infos.bool_tbl) {
        plot_bool(variant_infos.bool_tbl[k], k, variant_infos.violations, idxs);
    }
    let sample_bounds = get_sample_bounds()

    for (k in trios[0].tbl) {
        for (var i = 0; i < trios.length; i++) {
            var vals = trios[i].tbl[k].slice(0);
            plot_field(vals, k, k, true, trios[i].violations, idxs);
        }
    }
}

function initialize() {
    var seen = {}
    for (i = 0; i < variant_infos.filters.length; i++) {
        seen[variant_infos.filters[i]] = true
    };
    var nseen = 0;
    var haspass = false

    set_stats()

    // build global filters
    for (k in seen) { nseen += 1; if (k == "PASS") { haspass = true; } }
    let labels = Object.keys(seen)
    labels = labels.sort()
    if (nseen > 0) {
        jQuery('#filter-column').append(`<div class="btn-group-toggle" id="global-filters" data-toggle="buttons">`)
        labels.forEach(label => {
            let checked = ((!haspass) || (label == "PASS") ? "checked" : "");
            let active = checked == "checked" ? "active" : ""
            jQuery('#global-filters').append(`
                <label class="btn btn-outline-primary m-1 p-1 ${active}">
                    <input type="checkbox" autocomplete="off" class="slivar-filters" ${checked} name=${label} value=${label} id=${label}-check>
                        ${label}
                    </input>
                </label>
            `)
        });
        jQuery('#filter-column').append('</div>')
        jQuery('.slivar-filters').on('change', function () {
            jQuery('.slivar-checkbox').first().trigger('change');
        });
    }

    trios.forEach(function (trio) {
        sort_trio(trio)
    })

    // variants lengths filter
    add_slider(variant_infos.variant_lengths, "variant-length", "Variant Lengths", false, new Set(), variant_infos.violations)
    // variant filters: float
    for (k in variant_infos.float_tbl) {
        var vals = variant_infos.float_tbl[k];
        add_slider(vals, k, k, false, new Set(), variant_infos.violations)
    }
    // variant filters: bool
    for (k in variant_infos.bool_tbl) {
        let vals = variant_infos.bool_tbl[k];
        add_bool(vals, k, new Set(), variant_infos.violations)
    }

    for (k in trios[0].tbl) {
        add_slider(trios[0].tbl[k], k, k, true, new Set(), trios[0].violations)
        for (var i = 1; i < trios.length; i++) {
            var vals = trios[i].tbl[k];
            plot_field(vals, k, k, true, trios[i].violations, new Set());
        }
    }

    // whenever a checkbox (and slider) changes, re-draw all the plots
    // wrap in another function because update_info_plots expects first arg to be a
    // set. if empty, it creates it.
    jQuery(`.slivar-changer`).on('change', function (e, o) {
        // console.log("slivar-changer change", this, e, o, e.target);
        let idxs = new Set(get_passing_idxs())
        main_plot(idxs)
        update_info_plots(idxs)
    })
    // invert y-values in attribute ROC plots
    jQuery(`.slivar-inverter`).on('change', function (e, o) {
        let idxs = new Set(get_passing_idxs())
        plot_field(variant_infos.float_tbl[this.name], this.name, this.name, false, variant_infos.violations, idxs);
    })

    main_plot()
    update_stats([])
}

$(document).ready(function () {
    initialize()
    // refresh scrollspy
    $('body').scrollspy({
        target: '#side-filter-nav',
        offset: 515,
    })
})
