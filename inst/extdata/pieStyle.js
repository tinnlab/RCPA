// Settings

pieColors = [
    "#316b9d",
    "#fce397",
    "#99cc83",
    "#f77a65",
    "#a6a1d0",
    "#fea9c4",
    "#74e7bc",
    "#febb73",
    "#1db4db",
    "#ffc5a6",
    "#b6c9fa",
    "#ee5437"
]

notSignificantColor = "#eee"

// available data in node.data()
// ID: node ID
// label: node label
// size: node size
// nDE1: number of DE genes in result 1
// nDE2: number of DE genes in result 2
// ...
// isSig1: is significant in result 1
// isSig2: is significant in result 2
// ...
// pvalue1: pvalue in result 1
// pvalue2: pvalue in result 2
// ...
// nResult: number of results

vizmap = [
    {
        selector: "node",
        css: {
            label: 'data(label)',
            width: function (node) {
                return Math.pow(node.data("size"), 3 / 4)
            },
            height: function (node) {
                return Math.pow(node.data("size"), 3 / 4)
            },
            "background-color": "#ccc",
            "text-valign": "bottom",
            "text-halign": "center",
            "text-margin-y": "10px",
            "color": "#000",
            "font-size": "12px",
            "border-width": function (node) {
                let nDE = 0;
                for (let i = 1; i <= nResult; ++i) {
                    nDE += node.data('nDE' + i)
                }
                return Math.pow(nDE, 0.5) / nResult / 5
            },
            'border-color': "#eb4034",
            'pie-size': '100%'
        }
    },
    {
        selector: 'edge',
        css: {
            "line-color": "rgb(200, 200, 200)",
            "target-arrow-shape": "none",
            "target-arrow-color": "rgb(0, 0, 0)",
            "width": function (edge) {
                return Math.pow(edge.data("weight"), 0.5) / 2
            },
            'curve-style': 'bezier'
        }
    }
];

// End of settings

nResult = rcy.cy.nodes()[0].data("nResult");
pie = {};

for (let i = 1; i <= nResult; ++i) {
    let pieColor = pieColors[i - 1];

    pie[`pie-${i}-background-color`] = (node) => {
        let data = node.data();
        let isSig = data['isSig' + i];

        let color = isSig ? pieColor : notSignificantColor

        return color
    };
}

for (let i = 1; i <= nResult; ++i) {
    pie[`pie-${i}-background-size`] = 100 / nResult
}

vizmap[0].css = Object.assign(vizmap[0].css, pie);
