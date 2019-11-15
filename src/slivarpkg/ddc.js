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

function roc(values, violations, invert, abs) {
    var A = values.map(function(val, i) {
        if(invert) {val = -val};
        if(abs){val = Math.abs(val)};
        return [val, violations[i]]
    })
    A.sort(function(a, b) {
        return a[0] - b[0]
    })
    result = {x:[], y:[]}
    var tps = 0
    var fps = 0
    let N = A.length
    A.forEach(function(pair, i) {
        // FP
        if(pair[1]) { fps += 1; }
        else {tps += 1; }
        // TODO: scale by A.length
        // but might be good to leave as numbers for now.
        result.x.push(fps)
        result.y.push(tps)
    })
    return result
}

function add_slider(values, name, label) {

    jQuery('#variant-sliders').append(`<hr>
<label for="${name}-slider"><h3>${label}</h3></label>
<div id=${name}-slider name=${name}-slider></div>
<div class="col-8" id=${name}-hist-plot></div>
<!--- Joe Q. How to make these plots side-by-side. -->
<div class="col-8" id=${name}-roc-plot></div>
    `)
    var vlmin = Math.min.apply(null, values);
    var vlmax = Math.max.apply(null, values);
    sliders[`${name}`] = noUiSlider.create(document.getElementById(`${name}-slider`), {
    tooltips: [wNumb({decimals: 3}), wNumb({decimals: 3})],
    start: [vlmin, vlmax],
    connect: true, 
    range: {min:vlmin, max:vlmax},

})
    var bin_size = (vlmax - vlmin) / 100;
    var inh_vals = []
    var vio_vals = []
    var vios = variant_infos.violations.slice();
    for(var k=0; k < values.length; k++){
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

    var roc_tr = roc(values, vios, false, false);

    var traces = [
        {type: 'histogram', x: inh_vals, name:"transmitted", autobinx: false, xbins: {size: bin_size}, histnorm: "count"},
        {type: 'histogram', x: vio_vals, name:"violations", autobinx: false,  xbins: {size: bin_size}, histnorm: "count"},
    ];

    var roc_layout = {xaxis: {title:"mendelian violations"}, yaxis: {title:"transmitted variants"}}

    plots.hists[`${name}`] = Plotly.newPlot(`${name}-hist-plot`, traces, layout);
    plots.rocs[`${name}`] = Plotly.newPlot(`${name}-roc-plot`, [roc_tr], roc_layout);
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
        jQuery('#filter-column').append(`<input type='checkbox' ${checked} name=${k} value=${k} id=${k}>${k}<br>`)
    }
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




}());

var colors = ['rgba(55,126,184,0.7)', 'rgba(228,26,28,0.7)', 'rgba(77,175,74,0.7)', 'rgba(152,78,163,0.7)', 'rgba(255,127,0,0.7)', 'rgba(166,86,40,0.7)', 'rgba(247,129,191,0.7)']

var layout = {
    autosize: true,
    margin: {t: 30, pad: 0},
    xaxis: {
        title: "Mendelian violations",
        showspikes: true,
        spikethickness: 1,
    hoverformat: '.4f',
    },
    yaxis: {
        title: "Transmitted variants",
        showspikes: true,
        spikethickness: 1,
    hoverformat: '.4f',
    },
    hovermode: 'closest',
    showlegend: true,
    legend: {
        xanchor: "right",
        yanchor: "bottom",
        y: 0.1,
        x: 1,
        orientation: "h",
        borderwidth: 1,
        bordercolor: '#eeeeee'
    },
}

var layout = {
    autosize: true,
    margin: {t: 30, pad: 0},
    xaxis: {
        title: "Variant size-class",
        hoverformat: '.3f',
    },
    yaxis: {
        title: "False discovery rate",
        hoverformat: '.3f',
    },
    hovermode: 'closest',
    showlegend: true,
    legend: {
        x: 0.4,
        y: 0.9,
        borderwidth: 1,
        bordercolor: '#eeeeee'
    },

}



