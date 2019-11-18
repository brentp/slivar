var sliders = {samples:{}};
var plots = {hists: {}, rocs: {}}

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

function roc(values, violations, invert, abs, name, filters_to_keep) {
    let filters = variant_infos.filters
    var A = values.map(function(val, i) {
        if(invert) {val = -val};
        if(abs){val = Math.abs(val)};
        return [val, violations[i]]
    }).filter(function(val, i) {
        return filters_to_keep.includes(filters[i])
    })
    A.sort(function(a, b) {
        return a[0] - b[0]
    })
    result = {x:[0], y:[0], text: [""]}
    var tps = 0
    var fps = 0
    let N = A.length
    if(N == 0) {
        return result
    }
    var txt = (abs ? "Math.abs(": "") + "INFO." + name + (abs ? ") ": " ") + (invert ? ">" : "<")
    // we can't just iterate over A, we also have to track the change in
    // cutoff so we don't draw a smooth curve when the cutoff value hasn't
    // changed.
    var last_val = A[0][0] - 1;
    A.forEach(function(pair, i) {
        // FP
        if(pair[1]) { fps += 1; }
        else {tps += 1; }
        var val = pair[0]
        if (val == last_val) { return; }
        last_val = val;

        // TODO: scale by A.length
        // but might be good to leave as numbers for now.
        result.x.push(fps)
        result.y.push(tps)
        result.text.push(`${txt} ${val} gives FDR: ${(fps/(tps + fps)).toFixed(3)}`)
    })
    result.x.push(fps)
    result.y.push(tps)
    result.text.push(`${txt} ${A[A.length-1][0]} gives FDR: ${(fps/(tps + fps)).toFixed(3)}`)


    return result
}

function add_slider(values, name, label) {

    jQuery('#variant-sliders').append(`<hr>
<label for=${name}-flip>flip values</label>
<input type="checkbox" class=slivar-checkbox id=${name}-flip name=${name}-flip>
<label for=${name}-abs>take absolute value</label>
<input type="checkbox" class=slivar-checkbox id=${name}-abs name=${name}-abs>

<div id=${name}-slider name=${name}-slider></div>
<label for="${name}-slider"><h3>${label}</h3></label>
<div class="row pt-2 pb-2 bg-light">
<div class="col-6" id=${name}-hist-plot></div>
<div class="col-6" id=${name}-roc-plot></div>
</div>
    `)
    var vlmin = Math.min.apply(null, values);
    var vlmax = Math.max.apply(null, values);
    sliders[`${name}`] = noUiSlider.create(document.getElementById(`${name}-slider`), {
        tooltips: [wNumb({decimals: 3}), wNumb({decimals: 3})],
        start: [vlmin, vlmax],
        connect: true, 
        range: {min:vlmin, max:vlmax},
    })
    plot_info(values, name, label)

    jQuery(`#${name}-abs, #${name}-flip`).on('change', function() {
        plot_info(values, name, label);
    });

}

function get_passing_info_idxs() {
    var info_ranges = {}
    var info_bools = {} // TODO
    for(k in sliders) {
        if(k == "samples") { continue; }
        var tmp = sliders[k].get()
        info_ranges[k] = [parseInt(tmp[0]), parseInt(tmp[1])]
    }
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
        if(!skip) { idxs.push(i) }
    }
    return idxs
}

function main_plot() {
    // main-roc-plot
    console.time('main-plot')
    var traces = [];
    // we only plot points that pass the INFO filters
    let idxs = new Set(get_passing_info_idxs())
    console.log(idxs.length)
    trios.forEach(function(trio) {
        var trace = {x:[], y:[], text:[], name: trio.sample_id}
        for(var i=0; i < trio.variant_idxs.length; i++) {
            if(!idxs.has(trio.variant_idxs[i])){ continue; }
            trace.x.push(1)
            trace.y.push(1)


        }
        traces.push(trace)

    })
    console.timeEnd('main-plot')
}

function plot_info(values, name, label) {
    var vlmin = Math.min.apply(null, values);
    var vlmax = Math.max.apply(null, values);
    var filters_to_keep = jQuery('.slivar-filters:checked').map(function(v) { return this.name }).toArray();

    var bin_size = (vlmax - vlmin) / 100;
    var inh_vals = []
    var vio_vals = []
    var vios = variant_infos.violations.slice();
    let filters = variant_infos.filters;
    for(var k=0; k < values.length; k++){
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
    var do_flip = jQuery(`#${name}-flip`).is(":checked")
    // TODO: remove absolute value. not very helpful.
    var do_abs = jQuery(`#${name}-abs`).is(":checked")
    var roc_tr = roc(values, vios, do_flip, do_abs, name, filters_to_keep);

    var traces = [
        {type: 'histogram', x: inh_vals, name:"transmitted", autobinx: false, xbins: {size: bin_size}, histnorm: "count"},
        {type: 'histogram', x: vio_vals, name:"violations", autobinx: false,  xbins: {size: bin_size}, histnorm: "count"},
    ];

    var roc_layout = {xaxis: {title:"mendelian violations"}, yaxis: {title:"transmitted variants"}, 
        height:150, margin: {t: 0, b:20, r:0},
        legend: {x: 0.9, y: 0.9},
    }

    plots.hists[`${name}`] = Plotly.newPlot(`${name}-hist-plot`, traces, layout, {displayModeBar: false});
    plots.rocs[`${name}`] = Plotly.newPlot(`${name}-roc-plot`, [roc_tr], roc_layout, {displayModeBar: false});
}

function add_sample_slider(values, name) {
    jQuery('#sample-sliders').append(`
<label for="${name}-slider">${name}</label>
<div id=${name}-slider name=${name}-slider></div>
    `)
    var vlmin = Math.min.apply(null, values);
    var vlmax = Math.max.apply(null, values);
    sliders.samples[`${name}-slider`] = noUiSlider.create(document.getElementById(`${name}-slider`), {
        tooltips: [wNumb({decimals: 3}), wNumb({decimals: 3})],
    start: [vlmin, vlmax],
    connect: true, 
    range: {min:vlmin, max:vlmax},
})

}

(function(){
var seen = {}
for(i=0; i < variant_infos.filters.length;i++) {
    seen[variant_infos.filters[i]] = true
};
var nseen = 0;
var haspass = false
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

add_slider(variant_infos.variant_lengths, "variant-length", "variant lengths")

for(k in variant_infos.float_tbl){
    var vals = variant_infos.float_tbl[k];
    add_slider(vals, k, k)
}
// TODO: handle booleans

for(k in trios[0].kid_tbl) {
    var vals = trios[0].kid_tbl[k];
    add_sample_slider(vals, k)
}


main_plot()

}());

