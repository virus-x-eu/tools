$(document).ready(function () {
  const client = new $.es.Client({
    hosts: 'localhost'
  });
  const default_page_size = 100;

  $('#loading-indicator').hide();

  $('#search-button').click(function () {
    run_search();
  });
  $('#scroll-button').click(function () {
    if ((window.scroll_page * default_page_size) < window.total) {
      scroll_page();
    }
  });
  $('#export-button').click(function () {
    export_submission_ids();
  });
  $('#query-input').keypress(function (e) {
    if (e.which === 13) {
      run_search();
    }
  });

  function generate_search(user_search_text, page_size) {
    return {
      index: 'submission',
      body: {
        query: {
          multi_match: {
            query: user_search_text,
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
          },
          layout: {
            terms: {field: "experiment1.layout.keyword"}
          },
          library_source: {
            terms: {field: "experiment1.library_source.keyword"}
          }
        }
      },
      size: page_size,
      _source: ['sample_count', 'run_count'],
      scroll: '1h',
    };
  }

  async function run_search() {
    $('#loading-indicator').show();
    const response = await client.search(generate_search($('#query-input').val(), default_page_size));
    $('#loading-indicator').hide();
    window.scroll_id = response._scroll_id;
    window.scroll_page = 1;
    window.total = response.hits.total.value;

    $('#total-hits').text('' + response.hits.total.value + ' hits in ' + response.took + 'ms');
    $('#result-count').text('(' + response.hits.total.value + ')');

    update_table(response);

    for (k of ['platform', 'strategy', 'layout', 'library_source']) {
      $('#' + k + '-distribution').empty();
      for (const kx of response.aggregations[k].buckets) {
        $('#' + k + '-distribution').append('<li>' + kx.key + ' <b>(' + kx.doc_count + ')</b></li>');
      }
    }

    $('#export-count').text(response.hits.total.value);
  }

  async function scroll_page() {
    $('#loading-indicator').show();
    const response = await client.scroll({
      scrollId: window.scroll_id,
      scroll: '1h'
    });
    $('#loading-indicator').hide();
    window.scroll_page += 1;
    update_table(response);
  }

  function update_table(response) {
    $('#page-indicator').empty();
    $('#page-indicator').append('' + ((window.scroll_page - 1) * default_page_size + 1) + ' to ' + ((response.hits.total.value < default_page_size || ((window.scroll_page * default_page_size) > response.hits.total.value)) ? response.hits.total.value : (window.scroll_page * default_page_size)) + ' of ' + response.hits.total.value);
    $('#submissions-table').empty();
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
  }

  async function export_submission_ids() {
    $('#loading-indicator').show();
    const response = await client.search(generate_search($('#query-input').val(), 200));
    const all_submission_ids = [];
    const response_queue = [];

    response_queue.push(response);

    while (response_queue.length) {
      const body = response_queue.shift();

      for (const hit of body.hits.hits) {
        all_submission_ids.push(hit._id);
      }

      if (body.hits.total.value === all_submission_ids.length) {
        $('#loading-indicator').hide();
        break
      }

      response_queue.push(
        await client.scroll({
          scrollId: body._scroll_id,
          scroll: '30s',
        })
      );
    }

    var blob = new Blob([all_submission_ids.join('\n')], {type: 'text/plain'});
    $('#download-submission-ids').attr('href', (window.URL || window.webkitURL).createObjectURL(blob));
    $('#download-submission-ids')[0].click();
  }
});