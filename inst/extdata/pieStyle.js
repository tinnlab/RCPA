vizmap = [
    {
        selector: "node",
        css: {
            label: 'data(label)',
            width: 'data(size)',
            height: 'data(size)',
            "background-color": "#ccc",
            "text-valign": "bottom",
            "text-halign": "center",
            "text-margin-y": "3px",
            "color": "#000",
            "font-size": "12px",
            "border-width": 'data(borderWidth)',
            'border-color': 'data(borderColor)',
            'pie-size': '100%'
        }
    },
    {
        selector: 'edge',
        css: {
            "line-color": "data(edgeColor)",
            "target-arrow-shape": "none",
            "target-arrow-color": "rgb(0, 0, 0)",
            "width": 'data(weight)',
            'curve-style': 'bezier'
        }
    }
];

nResult = rcy.cy.nodes()[0].data("nResult");
pie = {};

for (let i = 1; i <= nResult; ++i) {
    pie[`pie-${i}-background-color`] = (node) => node.data("color" + i);
    pie[`pie-${i}-background-size`] = 100 / nResult
}

vizmap[0].css = Object.assign(vizmap[0].css, pie);
