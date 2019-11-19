var sliders = {samples:{}};
var plots = {hists: {}, rocs: {}}

var stats = {}

function set_stats() {
    stats.n_variants = variant_infos.violations.length;
    stats.n_violations = 0
    for(i=0;i<variant_infos.violations.length;i++){ stats.n_violations += variant_infos.violations[i]; }
}

function update_stats(idxs) {
    if(stats.filtered == undefined){ stats.filtered = {} }
    stats.filtered.n_violations = 0
    stats.filtered.n_variants = idxs.size
    idxs.forEach(function(e) {
        stats.filtered.n_violations += variant_infos.violations[e]
    })

    jQuery('#message').html(`<pre>
original variants: ${stats.n_variants}
variants after filtering: ${stats.filtered.n_variants} (${(100 * stats.filtered.n_variants / stats.n_variants ).toFixed(1)})

original violations: ${stats.n_violations}
violations after filtering: ${stats.filtered.n_violations} (${(100 * stats.filtered.n_violations / stats.n_violations ).toFixed(3)})

    </pre>`)
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
        to: function(a) {
            return a.toFixed(opts.decimals).replace(/([0-9]+(\.[0-9]+[1-9])?)(\.?0+$)/,'$1')
        },
        from: function(txt) {
            return parseFloat(txt);
        }
    };
}

var colors = ['rgba(55,126,184,0.7)', 'rgba(228,26,28,0.7)', 'rgba(77,175,74,0.7)', 'rgba(152,78,163,0.7)', 'rgba(255,127,0,0.7)', 'rgba(166,86,40,0.7)', 'rgba(247,129,191,0.7)']

function roc(values, violations, invert, name, filters_to_keep, idxs) {
    let filters = variant_infos.filters

    var Aorig = values.map(function(val, i) {
        if(invert) {val = -val};
        return [val, violations[i]]
    })
    var Afilt = Aorig.filter(function(val, i) {
        return (idxs.size == 0 || idxs.has(i)) && filters_to_keep.includes(filters[i])
    })
    Aorig.sort(function(a, b) {
        return a[0] - b[0]
    })
    Afilt.sort(function(a, b) {
        return a[0] - b[0]
    })
    let N = Aorig.length
    var txt = "INFO." + name + (invert ? ">" : "<")
    // we can't just iterate over A, we also have to track the change in
    // cutoff so we don't draw a smooth curve when the cutoff value hasn't
    // changed.
    result = [{x:[0], y:[0], text: [""], name: "unfiltered", marker: {color: colors[2]}},
              {x:[0], y:[0], text: [""], name: "filtered", marker: {color: colors[3]}}]
    if(N == 0) {
        return result
    }

    [Aorig, Afilt].forEach(function(A, Ai) {
        var last_val = A[0][0] - 1;
        var tps = 0
        var fps = 0
        A.forEach(function(pair, i) {
            // FP
            if(pair[1]) { fps += 1; }
            else {tps += 1; }
            var val = pair[0]
            if (val == last_val) { return; }
            last_val = val;

            // TODO: scale by A.length
            // but might be good to leave as numbers for now.
            result[Ai].x.push(fps)
            result[Ai].y.push(tps)
            result[Ai].text.push(`${txt} ${val} gives FDR: ${(fps/(tps + fps)).toFixed(3)}`)
        })
        result[Ai].x.push(fps)
        result[Ai].y.push(tps)
        result[Ai].text.push(`${txt} ${A[A.length-1][0]} gives FDR: ${(fps/(tps + fps)).toFixed(3)}`)
    })

    return result
}

function arr_min_max(arr) {
    var min = arr[0];
    var max = arr[0]
    for(i=1;i < arr.length; i++){
        if(arr[i] < min) {min = arr[i]; }
        if(arr[i] > max) {max = arr[i]; }
    }
    if(Math.abs(parseInt(min) - min) > 0.001) {
        min -= 0.01;
        max += 0.01;
    }
    return [min, max]
}


function add_bool(values, name, idxs){
    jQuery('#variant-sliders').append(`<hr>
<h3>${name}</h3>
<div class="row pt-2 pb-2 bg-light">
<label for=${name}-yes>yes</label>
<input type="checkbox" class="slivar-checkbox slivar-bool-checkbox slivar-changer" checked id=${name}-yes name=${name}-yes>
<label for=${name}-yes>no</label>
<input type="checkbox" class="slivar-checkbox slivar-bool-checkbox slivar-changer" checked id=${name}-no name=${name}-no>
</div>

<div class="row pt-2 pb-2 bg-light">
<div class="col-6" id=${name}-bar-plot></div>
</div>
    `)
    plot_bool(values, name, variant_infos.violations, idxs)
}

function add_slider(values, name, label, is_fmt_field, idxs) {

    var prefix = is_fmt_field ? "sample-": ""

    jQuery(is_fmt_field ? '#sample-sliders' : '#variant-sliders').append(`<hr>
<label for=${prefix}${name}-flip>flip values</label>
<input type="checkbox" class="slivar-checkbox slivar-changer" id=${prefix}${name}-flip name=${prefix}${name}-flip>

<div id=${prefix}${name}-slider name=${prefix}${name}-slider></div>
<label for="${prefix}${name}-slider"><h3>${label}</h3></label>
<div class="row pt-2 pb-2 bg-light">
<div class="col-6" id=${prefix}${name}-hist-plot></div>
<div class="col-6" id=${prefix}${name}-roc-plot></div>
</div>
    `)
    var [vlmin, vlmax] = arr_min_max(values); // can't use apply or destructring as we run out of stack
    var fmtF = name == "variant-length" ? iNumb : wNumb
    
    
    sliders[`${prefix}${name}`] = noUiSlider.create(document.getElementById(`${prefix}${name}-slider`), {
        tooltips: [fmtF({decimals: 3}), fmtF({decimals: 3})],
        start: [vlmin, vlmax],
        connect: true, 
        range: {min:vlmin, max:vlmax},
    })

    sliders[`${prefix}${name}`].on('change', function() {
        let idxs = new Set(get_passing_info_idxs())
        main_plot(idxs)
        update_info_plots(idxs)
    })

    plot_field(values, name, label, prefix, variant_infos.violations, idxs)
}

function get_sample_bounds() {
    var fmt_ranges = {}
    for(k in sliders) {
        if(k == "samples") { continue; }
        if(!k.startsWith("sample-")) {
            continue
        }
        let tmp = sliders[k].get()
        // remove sample- prefix
        fmt_ranges[k.substring(7)] = [parseFloat(tmp[0]), parseFloat(tmp[1])]
    }
    return fmt_ranges
}

function get_passing_info_idxs() {
    var info_ranges = {}
    var fmt_ranges = {}
    var info_bools = {} // TODO
    for(k in sliders) {
        if(k == "samples") { continue; }
        if(k.startsWith("sample-")) {
            continue
        } else {
            var tmp = sliders[k].get()
            info_ranges[k] = [parseFloat(tmp[0]), parseFloat(tmp[1])]
        }
    }

    jQuery('.slivar-bool-checkbox').each(function() {
        var n = this.name
        var i = -1;
        if(n.endsWith("-yes")){
            n = n.substring(0, n.length - 4)
            i = 1
        } else if(n.endsWith("-no")) {
            n = n.substring(0, n.length - 3)
            i = 0
        } else {
            alert("unknown bool type")
        }
        if(!(n in info_bools)) { info_bools[n] = [false, false]; }
            
        // so now, if info_bools[$x][0] is true, then we include things without
        // the flag. if [1] is true, we include things with the flag.
        info_bools[n][i] = this.checked

    })
    console.log(info_bools)


    var filters_to_keep = jQuery('.slivar-filters:checked').map(function(v) { return this.name }).toArray();

    var variant_lengths_ranges = info_ranges['variant-length']
    delete info_ranges['variant-length'];

    var idxs = []
    for(var i = 0; i < variant_infos.variant_lengths.length; i++){
        var vl = variant_infos.variant_lengths[i];
        if(vl < variant_lengths_ranges[0] || vl > variant_lengths_ranges[1]) { continue;}
        var skip = false;
        for(k in info_ranges){
            var vf = variant_infos.float_tbl[k][i];
            if (vf < info_ranges[k][0] || vf > info_ranges[k][1]){
                skip = true;
                break;
            }
        }
        if(skip) { continue}
        for(k in info_bools) {
            var fl = variant_infos.bool_tbl[k][i];
            if(!info_bools[k][Number(fl)]){
                skip = true;
                break;
            }
        }
        if(!skip) { idxs.push(i) }
    }
    return idxs
}

const ROC_VAR = "GQ"

function adjust_with_order(arr, order) {
    var result = new Array(arr.length)

    for(var i = 0; i < arr.length; i++){
        result[i] = arr[order[i]];
    }
    return result;
}

function sort_trio(trio) {
    // on page-load, we sort each trio so that the kid's lowest ABs
    // are first so we can draw the ROC curve
    var order = new Array(trio.violations.length)
    for(var i = 0; i < trio.violations.length; i++){ order[i] = i; }

    var roc_var = trio.tbl[ROC_VAR];
    //order.sort(function(a, b) { return roc_var[a] - roc_var[b]; });
    order.sort(function(a, b) { return (roc_var[b] - roc_var[a]); });

    for(k in trio.tbl) {
        trio.tbl[k] = adjust_with_order(trio.tbl[k], order)
    }
    trio.violations = adjust_with_order(trio.violations, order);
    trio.variant_idxs = adjust_with_order(trio.variant_idxs, order);
}


function sample_skippable(trio, sample_bounds, i) {
    for(k in sample_bounds) {
        let val = trio.tbl[k][i]
        if(val < sample_bounds[k][0] || val > sample_bounds[k][1]) { return true; }
    }
    return false;
}

function main_plot(idxs) {
    // main-roc-plot
    console.time('main-plot')
    var traces = [];
    // we only plot points that pass the INFO filters
    if(idxs == undefined) {
        idxs = new Set(get_passing_info_idxs())
    }
    let sample_bounds = get_sample_bounds()
    console.log(idxs.length)
    trios.forEach(function(trio) {
        var tps = 0
        var fps = 0
        var last_val = trio.tbl[ROC_VAR] - 1;
        var trace = {x:[], y:[], text:[], name: trio.sample_id}

        for(var i=0; i < trio.variant_idxs.length; i++) {
            if(!idxs.has(trio.variant_idxs[i])){ continue; }

            if(sample_skippable(trio, sample_bounds, i)){ continue; }
            var vio = trio.violations[i] //== 0 && trio.mom_alts[i] == 0 && trio.kid_alts[i] == 1;
            if(vio) { fps += 1; }
            else { tps += 1; }
            if(trio.tbl[ROC_VAR][i] == last_val) { continue; }
            last_val = trio.tbl[ROC_VAR][i];

            trace.x.push(fps);
            trace.y.push(tps);
            trace.text.push(`${ROC_VAR} > ${last_val} fpr: ${(fps / (fps + tps)).toFixed(3)}`)

        }
        console.log(`tps: ${tps} fps: ${fps} N: ${trio.variant_idxs.length}`)
        trace.x.push(fps);
        trace.y.push(tps);
        trace.text.push(`${ROC_VAR} > ${last_val} fpr: ${(fps / (fps + tps)).toFixed(3)}`)
        traces.push(trace)
    })
    console.timeEnd('main-plot')

    var layout = {xaxis: {title:"Violations"}, yaxis: {title: "Transmissions", autorange:true}}
    Plotly.newPlot('main-roc-plot', traces, layout);
}

function plot_bool(values, name, vios, idxs) {
    var filters_to_keep = jQuery('.slivar-filters:checked').map(function(v) { return this.name }).toArray();
    // convert bool to index (false==0) and increment counter for passing
    // variants
    var inh_counts = [0, 0]
    var vio_counts = [0, 0]
    
    var filt_inh_counts = [0, 0]
    var filt_vio_counts = [0, 0]
    let filters = variant_infos.filters;
    for(var k=0; k < values.length; k++){
        if(vios[k]) {
            vio_counts[Number(values[k])]++
        } else {
            inh_counts[Number(values[k])]++
        }
        if(idxs.size > 0 && !idxs.has(k)) { continue; }
        if(!filters_to_keep.includes(filters[k])) { continue; }
        if(vios[k]) {
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
    if(Y) {
        vio[0] += vio_counts[1]
        tra[0] += inh_counts[1]

        vio[1] += filt_vio_counts[1]
        tra[1] += filt_inh_counts[1]
    }
    if(N) {
        vio[0] += vio_counts[0]
        tra[0] += inh_counts[0]

        vio[1] += filt_vio_counts[0]
        tra[1] += filt_inh_counts[0]
    }

    let xlabels = ["unfiltered", "filtered"]
    var traces = [
                  {y: vio, x: xlabels, name:"violations", type: "bar"},
                  {y: tra, x: xlabels, name:"transmissions", type: "bar"},
                ]
    var layout = {barmode: 'group'};
    Plotly.newPlot(`${name}-bar-plot`, traces, layout)


}

function plot_field(values, name, label, prefix, vios, idxs) {
    var [vlmin, vlmax] = arr_min_max(values);
    var filters_to_keep = jQuery('.slivar-filters:checked').map(function(v) { return this.name }).toArray();

    var inh_vals = []
    var vio_vals = []
    let filters = variant_infos.filters;
    for(var k=0; k < values.length; k++){
        if(idxs.size > 0 && !idxs.has(k)) { continue; }
        if(!filters_to_keep.includes(filters[k])) { continue; }
        if(vios[k]) {
            vio_vals.push(values[k])
        } else {
            inh_vals.push(values[k])
        }
    }

    var layout = {barmode:"overlay", xaxis: {}, yaxis: {title: "Count", autorange:true, type:"log"}, 
        height:150, margin: {t: 0, b:20, r:0},
        legend: {x: 0.9, y: 0.9},
         
    };
    var do_flip = jQuery(`#${prefix}${name}-flip`).is(":checked")
    var roc_trs = roc(values, vios, do_flip, name, filters_to_keep, idxs);
    var bin_size = (vlmax - vlmin) / 100;

    var traces = [
        {type: 'histogram', x: inh_vals, name:"transmitted", autobinx: false, xbins: {size: bin_size}, histnorm: "count"},
        {type: 'histogram', x: vio_vals, name:"violations", autobinx: false,  xbins: {size: bin_size}, histnorm: "count"},
    ];

    var roc_layout = {xaxis: {title:"mendelian violations"}, yaxis: {title:"transmitted variants"}, 
        height:150, margin: {t: 0, b:20, r:0},
        legend: {x: 0.9, y: 0.9},
    }

    plots.hists[`${prefix}${name}`] = Plotly.newPlot(`${prefix}${name}-hist-plot`, traces, layout, {displayModeBar: false});
    plots.rocs[`${prefix}${name}`] = Plotly.newPlot(`${prefix}${name}-roc-plot`, roc_trs, roc_layout, {displayModeBar: false});
}

(function(){
var seen = {}
for(i=0; i < variant_infos.filters.length;i++) {
    seen[variant_infos.filters[i]] = true
};
var nseen = 0;
var haspass = false

set_stats()

for(k in seen) { nseen += 1; if(k == "PASS"){ haspass = true; }}
if(nseen > 0) {
    for(k in seen){
    var checked = ((!haspass) || (k == "PASS") ? "checked": "");
        jQuery('#filter-column').append(`<input type='checkbox' class=slivar-filters ${checked} name=${k} value=${k} id=${k}-check >${k}<br>`)
    }
    jQuery('.slivar-filters').on('change', function() {
        jQuery('.slivar-checkbox').trigger('change');
    });
}

console.time('trio-sort')
trios.forEach(function(trio) {
    sort_trio(trio)
})
console.timeEnd('trio-sort')

add_slider(variant_infos.variant_lengths, "variant-length", "variant lengths", false, new Set())

for(k in variant_infos.float_tbl){
    var vals = variant_infos.float_tbl[k];
    add_slider(vals, k, k, false, new Set())
}

for(k in variant_infos.bool_tbl) {
    let vals = variant_infos.bool_tbl[k];
    add_bool(vals, k, new Set())
}

// TODO: handle booleans

for(k in trios[0].tbl) {
    var vals = trios[0].tbl[k];
    add_slider(vals, k, k, true, new Set())
}

// whenever a checkbox (and slider) changes, re-draw all the plots
// wrap in another function because update_info_plots expects first arg to be a
// set. if empty, it creates it.
jQuery(`.slivar-changer`).on('change', function() { update_info_plots() })


main_plot()

}());

function update_info_plots(idxs) {
    if(idxs == undefined) {
        idxs = new Set(get_passing_info_idxs())
    } else {
        console.log("got idxs", idxs)
    }
    let prefix = "";
    update_stats(idxs)
    plot_field(variant_infos.variant_lengths, "variant-length", "variant-lengths", prefix, variant_infos.violations, idxs);
    for(k in variant_infos.float_tbl) {
        plot_field(variant_infos.float_tbl[k], k, k, prefix, variant_infos.violations, idxs);
    }
    for(k in variant_infos.bool_tbl) {
        plot_bool(variant_infos.bool_tbl[k], k, variant_infos.violations, idxs);
    }
}
