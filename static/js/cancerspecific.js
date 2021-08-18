url = new URL(window.location.href)

$.ajax({
    url: "/CancerSpecific",
    type: "GET",
    contentType: 'application/json;charset=UTF-8',
    data: {
        cancer: url.searchParams.get("cancer"),
    },
    dataType: "json",
    success: function(data) {
        Plotly.newPlot('sample_graph', data, {responsive: true});
    }
});