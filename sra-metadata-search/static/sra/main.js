$(document).ready(function () {
  const client = new $.es.Client({
    hosts: 'localhost'
  });

  $('#loading-indicator').hide();

  $('#search-button').click(function () {
    run_search();
  });
  $('#query-input').keypress(function (e) {
    if (e.which === 13) {
      run_search();
    }
  });

  async function run_search() {
    $('#loading-indicator').show();
    const response = await client.search({
      index: 'submission',
      body: {
        query: {
          multi_match: {
            query: $('#query-input').val(),
            operator: "and",
            fields: [
              'sample1.*',
              'experiment1.*'
            ],
          }
        },
        highlight: {
          fields: {
            "*": {}
          }
        },
        aggs: {
          platform: {
            terms: {field: "experiment1.platform.keyword"}
          },
          strategy: {
            terms: {field: "experiment1.strategy.keyword"}
          }
        }
      },
      size: 1000,
      _source: ['sample_count', 'run_count']
    });
    $('#loading-indicator').hide();

    $('#total-hits').text('' + response.hits.total.value + ' hits in ' + response.took + 'ms');
    $('#result-count').text('(' + response.hits.total.value + ')');

    $('#submissions-table').empty();
    $('#platform-distribution').empty();
    $('#strategy-distribution').empty();

    for (const hit of response.hits.hits) {
      let highlights = "";
      if ('highlight' in hit) {
        for (const field of Object.keys(hit.highlight)) {
          for (const match of hit.highlight[field]) {
            highlights += '<div style="margin-bottom: 10px;"><span>' + match + '</span> <sup style="opacity: 0.5">' + field + '</sup></div>';
          }
        }
      }
      const st = $('#submissions-table');
      st.append('<tr>');
      st.append('<td style="font-family: monospace">' + hit._id + '</td>');
      st.append('<td>' + highlights + '</td>');
      st.append('<td style="font-family: monospace">' + (hit._source.sample_count ? hit._source.sample_count : '-') + '</td>');
      st.append('<td style="font-family: monospace">' + (hit._source.run_count ? hit._source.run_count : '-') + '</td>');
      st.append('</tr>');
    }

    for (const platform of response.aggregations.platform.buckets) {
      $('#platform-distribution').append('<li>' + platform.key + ' <b>(' + platform.doc_count + ')</b></li>');
    }
    for (const strategy of response.aggregations.strategy.buckets) {
      $('#strategy-distribution').append('<li>' + strategy.key + ' <b>(' + strategy.doc_count + ')</b></li>');
    }
  }
});