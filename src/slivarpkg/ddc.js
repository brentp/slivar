var ROC_VAR = "GQ"
var ROC_SIGN = ROC_VAR == "AB" ? "<" : ">"
var ranges = { samples: {} }
var plots = { hists: {}, rocs: {} }
var stats = {filtered: {}}
var sorted = {}
// track select/de-select events
var applied_filters = []
const colors = ["#4e79a7", "#f28e2c", "#e15759", "#76b7b2", "#59a14f", "#edc949", "#af7aa1", "#ff9da7", "#9c755f", "#bab0ab"]
const primary_color = "#2F3C48"
const main_roc_layout = {
    height: 400,
    margin: { t: 20, r: 0, b: 0, l: 60 },
    xaxis: { automargin: true, title: { text: "Mendelian Violations", standoff: 100 }, rangemode: 'tozero', },
    yaxis: { automargin: true, title: { text: "Transmitted Variants", standoff: 10 }, rangemode: 'tozero', },
    font: { size: 15 }
}
const aux_plot_height = 250
const histogram_layout = {
    barmode: "overlay",
    xaxis: { automargin: true, fixedrange: true }, //, rangeslider: {} },
    yaxis: { title: { text: "Count", standoff: 20 }, autorange: true, automargin: true, type: "log", fixedrange: true },
    height: aux_plot_height,
    margin: { t: 10, b: 40, r: 0, l: 50, pad: 0 },
    legend: {
        xanchor: "right",
        yanchor: "top",
        y: 1.1,
        x: 1,
        orientation: "h",
    },
    dragmode: 'select',
    selectdirection: "h"
}
const roc_layout = {
    xaxis: { automargin: true, title: { text: "Violations", standoff: 20 } },
    yaxis: { automargin: true, title: { text: "Transmitted", standoff: 20 }, rangemode: 'tozero' },
    hovermode: "y",
    hoverinfo:"text",
    height: aux_plot_height,
    margin: { t: 10, b: 20, r: 20, l: 30, pad: 0 },
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
    var filters_to_keep = jQuery('.slivar-filters:checked').map(function (v) { return this.name }).toArray();
    var filters = variant_infos.filters

    stats.n_violations = 0
    stats.n_transmitted = 0
    for (let i = 0; i < variant_infos.violations.length; i++) {
        if (!filters_to_keep.includes(filters[i])) { continue; }
        stats.n_violations += variant_infos.violations[i];
        stats.n_transmitted += (1 - Number(variant_infos.violations[i]));
    }
}

function update_stats(idxs) {
    let fpr = (stats.n_violations / (stats.n_violations + stats.n_transmitted)).toFixed(3)
    stats.filtered.n_violations = 0
    stats.filtered.n_variants = idxs.size
    stats.filtered.n_transmitted = 0
    idxs.forEach(function (e) {
        stats.filtered.n_violations += variant_infos.violations[e]
    })
    stats.filtered.n_transmitted = idxs.size - stats.filtered.n_violations
    // no filters have been applied
    jQuery('#original-variants').html(`${stats.n_transmitted}`)
    jQuery('#original-violations').html(`${stats.n_violations}`)
    jQuery('#original-fpr').html(fpr)
    if(idxs.size != 0){
        var f_fpr = 0
        if (stats.filtered.n_violations + stats.filtered.n_transmitted > 0) {
            f_fpr = (stats.filtered.n_violations / (stats.filtered.n_violations + stats.filtered.n_transmitted)).toFixed(3)
        }
        jQuery('#filtered-variants').html(`${stats.filtered.n_transmitted} <span class="text-info">(${(100 * stats.filtered.n_transmitted / stats.n_transmitted).toFixed(2)}%)</span>`)
        jQuery('#filtered-violations').html(`${stats.filtered.n_violations} <span class="text-info">(${(100 * stats.filtered.n_violations / stats.n_violations).toFixed(2)}%)</span>`)
        jQuery('#filtered-fpr').html(`${f_fpr}`)
    }

    jQuery('#filtering-on').attr('hidden', applied_filters.length == 0)
    if(applied_filters.length == 0){
                jQuery('#original-row-label').html('')
    }
    jQuery('#filtered-row').attr('hidden', applied_filters.length == 0)
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
            return [val, violations[i], i]
        })
        sorted[key].sort(function (a, b) {
            return a[0] - b[0]
        })
    }
    var Aorig = sorted[key];
    var Afilt = Aorig.filter(function (pair, pi) {
        let i = pair[2];

        return (idxs.size == 0 || idxs.has(i)) && filters_to_keep.includes(filters[i])
    })
    //console.timeEnd("deco-sort " + name)
    let N = Aorig.length
    var txt = "INFO." + name + " " + (invert ? ">" : "<") + " "
    if(name == "variant-length") {
        txt = name + " " + (invert ? ">" : "<") + " "
    }
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
                result[Ai].text.push(`${txt} ${val}<br>FDR: ${(fps / (tps + fps)).toFixed(3)}`)
            }
        })
        result[Ai].x.push(fps)
        result[Ai].y.push(tps)
        //console.log((fps / (tps + fps)).toFixed(3), "fps:", fps, " AI:", Ai)
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
            <div class="row">
                <div class="col-12 p-1" id=${name}-bar-plot></div>
            </div>
            <!-- Boolean selection -->
            <div class="row pb-3">
                <div class="col-1 pr-0">
                    <span class="selection-label">Include:</span>
                </div>
                <div class="col-11">
                    <div class="btn-group-toggle" data-toggle="buttons">
                        <label class="btn btn-outline-primary m-1 p-1 active cbtn">
                            <input type="checkbox" autocomplete="off" class="slivar-checkbox slivar-bool-checkbox slivar-changer" checked name=${name}-yes id=${name}-yes>
                                ${name}
                            </input>
                        </label>
                        <label class="btn btn-outline-primary m-1 p-1 active cbtn">
                            <input type="checkbox" autocomplete="off" class="slivar-checkbox slivar-bool-checkbox slivar-changer" checked name=${name}-no id=${name}-no>
                                not ${name}
                            </input>
                        </label>
                    </div>
                </div>
            </div>
        </div>
    `)
    plot_bool(values, name, variant_infos.violations)

    jQuery(`#${name}-yes,#${name}-no`).change(function () {
        plot_bool(values, name, variant_infos.violations)
    })
}

function add_slider(values, name, label, is_fmt_field, violations) {
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
                <div class="row">
                    <div class="col-12 p-1" id=${prefix}${name}-hist-plot></div>
                </div>
                <div class="row">
                    <div class="col-8 text-monospace" id="${prefix}${name}-code-field"></div>
                    <div class="col-4 text-right"></div>
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
                <div class="row">
                    <div class="col-8 p-1" id=${prefix}${name}-hist-plot></div>
                    <div class="col-4 p-1" id=${prefix}${name}-roc-plot></div>
                </div>
                <div class="row">
                    <div class="col-8 text-monospace" id="${prefix}${name}-code-field">
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

    plot_field(values, name, label, is_fmt_field, violations)
}

function get_sample_bounds() {
    var sample_bounds = {}
    for (k in ranges) {
        if (k == "samples") { continue; }
        if (!k.startsWith("sample-")) {
            continue
        }
        let tmp = ranges[k]
        // don't do anything if it's the same ranges as when the plot loaded.
        if (Math.abs(tmp.min - ranges[k].omin) < 0.005 && Math.abs(tmp.max - ranges[k].omax) < 0.005) {
            continue
        }
        // remove sample- prefix
        sample_bounds[k.substring(7)] = [tmp.min, tmp.max]
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
                    result.push(trio.variant_idxs[i])
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

    for (k in ranges) {
        if (k == "samples") { continue; }
        if (k.startsWith("sample-")) {
            continue
        } else {
            var tmp = ranges[k]
            if (Math.abs(tmp.min - ranges[k].omin) < 0.005 && Math.abs(tmp.max - ranges[k].omax) < 0.005) {
                continue
            }
            info_bounds[k] = [ranges[k].min, ranges[k].max]
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
    if(!(ROC_VAR in trio.tbl)){
        console.log("can't use {ROC_VAR} as it's not present")
        ROC_VAR = Object.keys(trio.tbl)[0]
        ROC_SIGN = ROC_VAR == "AB" ? "<" : ">"
    }

    var roc_var = trio.tbl[ROC_VAR];
    if (ROC_VAR == "AB") {
        var rv2 = new Array(roc_var.length)
        for (var i = 0; i < rv2.length; i++) {
            rv2[i] = 0.5 - Math.abs(0.5 - roc_var[i])
        }
        roc_var = rv2
    }
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
    if(idxs != undefined) {
        console.log("idxs:", idxs.size)
    }
    var traces = [];
    // we only plot points that pass the INFO filters
    if (idxs == undefined || idxs.size == 0) {
        idxs = new Set(get_passing_idxs())
    }
    update_stats(idxs)

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
            trace.text.push(`${ROC_VAR} ${ROC_SIGN} ${last_val} fpr: ${(fps / (fps + tps)).toFixed(3)}`)

        }
        trace.x.push(fps);
        trace.y.push(tps);
        trace.text.push(`${ROC_VAR} > ${last_val} fpr: ${(fps / (fps + tps)).toFixed(3)}`)
        traces.push(trace)
    })
    console.timeEnd('main-plot')

    Plotly.newPlot('main-roc-plot', traces, main_roc_layout)
}

function plot_bool(values, name, vios) {
    var filters_to_keep = jQuery('.slivar-filters:checked').map(function (v) { return this.name }).toArray();
    // convert bool to index (false==0) and increment counter for passing
    // variants
    var inh_counts = [0, 0]
    var vio_counts = [0, 0]

    var filt_inh_counts = [0, 0]
    var filt_vio_counts = [0, 0]
    let filters = variant_infos.filters;
    var Y = jQuery(`#${name}-yes`).is(":checked")
    var N = jQuery(`#${name}-no`).is(":checked")
    applied_filters = applied_filters.filter(i => ((i !== `${name}-yes`) && (i !== `${name}-no`)))
    // if both are selected it's same as not filtering
    if(Y && !N){
       applied_filters.push(`${name}-yes`)
    }
    if(N && !Y){
       applied_filters.push(`${name}-no`)
    }


    let info_bools = [N, Y];

    for (var k = 0; k < values.length; k++) {
        if (vios[k]) {
            vio_counts[Number(values[k])]++
        } else {
            inh_counts[Number(values[k])]++
        }
        if (!info_bools[Number(values[k])]) { continue; }
        //if (idxs.size > 0 && !idxs.has(k)) { continue; }
        if (!filters_to_keep.includes(filters[k])) { continue; }
        if (vios[k]) {
            filt_vio_counts[Number(values[k])]++
        } else {
            filt_inh_counts[Number(values[k])]++
        }
    }
    //console.log(vio_counts, filt_vio_counts)
    //console.log(Y, N)

    // unfilt, filt
    var vio = [0, 0]
    var tra = [0, 0]
    vio[0] += vio_counts[1]
    tra[0] += inh_counts[1]

    vio[1] += filt_vio_counts[1]
    tra[1] += filt_inh_counts[1]

    vio[0] += vio_counts[0]
    tra[0] += inh_counts[0]

    vio[1] += filt_vio_counts[0]
    tra[1] += filt_inh_counts[0]

    let xlabels = ["Unfiltered", "Filtered"]
    var traces = [
        { x: tra, y: xlabels, name: "Transmitted", type: "bar", orientation: 'h', marker: { color: colors[0] }, },
        { x: vio, y: xlabels, name: "Violations", type: "bar", orientation: 'h', marker: { color: colors[1] }, },
    ]
    Plotly.newPlot(`${name}-bar-plot`, traces, box_layout)
}

function plot_field(values, name, label, is_fmt_field, vios) {
    // update_only is used only when is_fmt_field to add samples to an existing
    // fmt roc plot
    var [vlmin, vlmax] = arr_min_max(values);
    //console.log("name:", name, [vlmin, vlmax])
    var filters_to_keep = jQuery('.slivar-filters:checked').map(function (v) { return this.name }).toArray();

    var inh_vals = []
    var vio_vals = []
    let filters = variant_infos.filters;
    //console.time("plot_field filtering:" + name)
    for (var k = 0; k < values.length; k++) {
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
        roc_trs = roc(values, vios, do_flip, name, filters_to_keep, new Set());
    }

    var traces = [
        { type: 'histogram', x: inh_vals, name: "Transmitted", autobinx: false, xbins: { size: bin_size }, histnorm: "count", marker: { color: colors[0] } }, //, opacity: 0.5 },
        { type: 'histogram', x: vio_vals, name: "Violations", autobinx: false, xbins: { size: bin_size }, histnorm: "count", marker: { color: colors[1] } }, //, opacity: 0.5 },
    ];

    let div = document.getElementById(`${prefix}${name}-hist-plot`)
    plots.hists[`${prefix}${name}`] = Plotly.newPlot(div, traces, histogram_layout, { displayModeBar: false });
    ranges[`${prefix}${name}`] = { min: vlmin, max: vlmax, omin: vlmin, omax: vlmax }

    let update = function (e) {
        if (this.id == `${prefix}${name}-btn-close`) {
            e = null
            jQuery(`#${prefix}${name}-btn-close`).remove()
            jQuery(`#${prefix}${name}-btn-lbl`).remove()
            // simulates plotly deselect event when clicking button
            Plotly.restyle(`${prefix}${name}-hist-plot`, {selectedpoints: [null]});
        }
        // single-click
        if (e === undefined) { return; }
        // de-select event; dbl-click to reset plot doesn't have range attr
        if (e === null || !("range" in e)) {
            // pop filter from array
            applied_filters = applied_filters.filter(i => i !== `${prefix}${name}`)
            // hide filtering labels if no filters are active
            if (applied_filters.length == 0) {
                jQuery('#filtering-on').attr('hidden', true)
                jQuery('#original-row-label').html('')
                jQuery('#filtered-row').attr('hidden', true)
            }
            // reset the ranges
            ranges[`${prefix}${name}`].min = ranges[`${prefix}${name}`].omin
            ranges[`${prefix}${name}`].max = ranges[`${prefix}${name}`].omax
        // (re-)select event
        } else {
            // change button label
            if(label != "Variant Lengths") {
                jQuery(`#${prefix}${name}-btn-lbl`).html(`${label} (${e.range.x[0].toFixed(1)} - ${e.range.x[1].toFixed(1)})`)
                // show select expression
                jQuery(`#${prefix}${name}-code-field`).html(`
                    <span class="small" id="${prefix}${name}-expression">INFO.${label} >= ${e.range.x[0].toFixed(3)} && INFO.${label} <= ${e.range.x[1].toFixed(3)}</span>
                    <button type="button" class="btn btn-sm" id="${prefix}${name}-btn-copy" data-name="${prefix}${name}" data-toggle="tooltip" data-placement="top" title="Copy to clipboard">
                        <i class="far fa-copy"></i>
                    </button>
                `)
                jQuery(`#${prefix}${name}-btn-copy`).on('click', copy_code_function)
                jQuery(`#${prefix}${name}-btn-copy`).tooltip('enable')
            } else {
                jQuery(`#${prefix}${name}-btn-lbl`).html(`${label} (${e.range.x[0].toFixed(1)} - ${e.range.x[1].toFixed(1)})`)
            }
            // show filtering label
            jQuery(`#filtering-on`).attr('hidden', false)
            jQuery('#original-row-label').html('Original')
            jQuery('#filtered-row').attr('hidden', false)
            applied_filters.push(`${prefix}${name}`)
            // show the button
            jQuery(`#${prefix}${name}-btn`).attr('hidden', false)
            ranges[`${prefix}${name}`].min = e.range.x[0]
            ranges[`${prefix}${name}`].max = e.range.x[1]
        }
        main_plot(new Set(get_passing_idxs()))
        if (!is_fmt_field) {
            var vmin = ranges[`${prefix}${name}`].min;
            var vmax = ranges[`${prefix}${name}`].max;

            var marginal_idxs = [];
            var filters_to_keep = jQuery('.slivar-filters:checked').map(function (v) { return this.name }).toArray();
            for (var k = 0; k < values.length; k++) {
                if (!filters_to_keep.includes(filters[k])) { continue; }
                let v = values[k]
                if (v >= vmin && v <= vmax) {
                    marginal_idxs.push(k);
                }
            }

            let do_flip = jQuery(`#${prefix}${name}-flip`).is(":checked")
            roc_trs = roc(values, vios, do_flip, name, filters_to_keep, new Set(marginal_idxs));
            plots.rocs[`${prefix}${name}`] = Plotly.newPlot(`${prefix}${name}-roc-plot`, roc_trs, roc_layout, { displayModeBar: false });
        }
    }

    // create filtering button for histogram
    jQuery('#active-filters').append(`
        <div class="btn-group mb-1 mr-1" id="${prefix}${name}-btn" hidden>
            <button type="button" class="btn btn-outline-primary btn-sm cbtn" id="${prefix}${name}-btn-lbl" onclick="window.location.href='#${prefix}${name}-wrapper';">
                ${label}
            </button>
            <button type="button" class="btn btn-primary btn-sm" id="${prefix}${name}-btn-close" data-name="${prefix}${name}" title="Remove filter">
                <i class="fas fa-times fa-lg"></i>
            </button>
        </div>
    `)
    // add click functions
    jQuery(`#${prefix}${name}-btn-close`).on('click', update)
    div.on('plotly_selected', update);
    div.on('plotly_deselect', update);
    // plot the roc curve next to histogram
    if (!is_fmt_field) {
        plots.rocs[`${prefix}${name}`] = Plotly.newPlot(`${prefix}${name}-roc-plot`, roc_trs, roc_layout, { displayModeBar: false });
    }
}

function copy_code_function() {
    let label = this.getAttribute('data-name')
    let $temp = $("<input>")
    $("body").append($temp)
    $temp.val($(`#${label}-expression`).text()).select()
    document.execCommand("copy")
    $temp.remove()
    // notify action
    jQuery(`#${label}-btn-copy`).attr("title", "Copied!").tooltip("_fixTitle").tooltip("show").attr("title", "Copy to clipboard").tooltip("_fixTitle");
}

function initialize_once() {
    var seen = {}
    for (i = 0; i < variant_infos.filters.length; i++) {
        seen[variant_infos.filters[i]] = true
    };
    var nseen = 0;
    var haspass = false


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
                <label class="btn btn-outline-primary m-1 p-1 ${active} cbtn">
                    <input type="checkbox" autocomplete="off" class="slivar-filters" ${checked} name=${label} value=${label} id=${label}-check>
                        ${label}
                    </input>
                </label>
            `)
        });
        jQuery('#filter-column').append('</div>')
    }
    jQuery('.slivar-filters').on('change', redraw)

    trios.forEach(function (trio) {
        sort_trio(trio)
    })
}

function redraw() {
    // variants lengths filter
    set_stats()
    jQuery('#sample-filter-nav').empty()
    jQuery('#sample-sliders').empty()
    jQuery('#variant-filter-nav').empty()
    jQuery('#variant-sliders').empty()

    // rebuild anchors and titles
    jQuery('#variant-sliders').append(`
        <a id="variant-filter-label"></a>
        <h1>Variant Filters</h1>
        (drag to select a region in each plot to choose cutoffs)
    `)
    jQuery('#sample-sliders').append(`
        <a id="sample-filter-label"></a>
        <h1>Sample Filters</h1>
    `)

    add_slider(variant_infos.variant_lengths, "variant-length", "Variant Lengths", false, variant_infos.violations)
    // variant filters: float
    for (k in variant_infos.float_tbl) {
        var vals = variant_infos.float_tbl[k];
        add_slider(vals, k, k, false, variant_infos.violations)
    }
    // variant filters: bool
    for (k in variant_infos.bool_tbl) {
        let vals = variant_infos.bool_tbl[k];
        add_bool(vals, k, new Set(), variant_infos.violations)
    }

    for (k in trios[0].tbl) {
        add_slider(trios[0].tbl[k], k, k, true, trios[0].violations)
        for (var i = 1; i < trios.length; i++) {
            var vals = trios[i].tbl[k];
            plot_field(vals, k, k, true, trios[i].violations);
        }
    }

    // whenever a checkbox (and slider) changes, re-draw all the plots
    // wrap in another function because update_info_plots expects first arg to be a
    // set. if empty, it creates it.
    jQuery(`.slivar-changer`).on('change', function (e, o) {
        // console.log("slivar-changer change", this, e, o, e.target);
        let idxs = new Set(get_passing_idxs())
        main_plot(idxs)
    })
    // invert y-values in attribute ROC plots
    jQuery(`.slivar-inverter`).on('change', function (e, o) {
        plot_field(variant_infos.float_tbl[this.name], this.name, this.name, false, variant_infos.violations);
    })

    main_plot()
}


$(document).ready(function () {
    initialize_once()
    redraw()
    // refresh scrollspy
    $('body').scrollspy({
        target: '#side-filter-nav',
        offset: 515,
    })
    // initialize tooltips
    $('[data-toggle="tooltip"]').tooltip()
})
