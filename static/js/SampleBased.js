url = new URL(window.location.href)

$.ajax({
    url: "/SampleBased",
    type: "GET",
    contentType: 'application/json;charset=UTF-8',
    data: {
        CanSampleId: url.searchParams.get("sampleid"),
        genename: url.searchParams.get("gene"),
        tissuetype: url.searchParams.get("tissue")
    },
    dataType: "json",
    success: function(data) {
        Plotly.newPlot('bargraph', data, {responsive: true});
    }
});
