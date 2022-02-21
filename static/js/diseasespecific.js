url = new URL(window.location.href)

$.ajax({
    url: "/DiseaseSpecific",
    type: "GET",
    contentType: 'application/json;charset=UTF-8',
    data: {
        disease: url.searchParams.get("disease"),
    },
    dataType: "json",
    success: function(data) {
        Plotly.newPlot('sample_graph', data, {responsive: true});
    }
});