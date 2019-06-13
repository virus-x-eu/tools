$(document).ready(function () {
  const client = new $.es.Client({
    hosts: 'localhost:' + window.location.port
  });
  const default_page_size = 50;

  const nf = new Intl.NumberFormat('en-US');

  $('#loading-indicator').hide();
  $('#query-input').focus();

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
    let query = {
      multi_match: {
        query: user_search_text,
        operator: "and",
        type: "cross_fields",
        fields: [
          'id',
          'sample1.*',
          'experiment1.*'
        ],
      }
    };
    if (user_search_text === "") {
      query = {match_all: {}}
    }
    let source = ['sample_count', 'run_count', 'experiment1_title', 'sample1_title', 'date'];
    if (page_size >= 200) {
      source = false;
    }
    return {
      index: 'submission',
      body: {
        query: query,
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
          },
          worldmap: {
            geotile_grid: {
              field: "sample1_location",
              precision: 8,
              size: 10000
            },
          },
          run_count: {
            sum: {
              field: "run_count"
            }
          },
          date: {
            date_histogram: {
              field: "date",
              interval: "year",
            }
          }
        }
      },
      size: page_size,
      _source: source,
      scroll: '20m',
    };
  }

  async function run_search() {
    $('#loading-indicator').show();
    const response = await client.search(generate_search($('#query-input').val(), default_page_size));
    $('#loading-indicator').hide();
    window.scroll_id = response._scroll_id;
    window.scroll_page = 1;
    window.total = response.hits.total.value;

    update_worldmap(response);

    $('#total-hits').empty();
    $('#total-hits').append('<b>' + nf.format(response.hits.total.value) + ' hits</b> in ' + response.took + 'ms');
    $('#result-count').text('(' + nf.format(response.hits.total.value) + ')');

    update_table(response);

    for (k of ['platform', 'strategy', 'layout', 'library_source']) {
      $('#' + k + '-distribution').empty();
      for (const kx of response.aggregations[k].buckets) {
        $('#' + k + '-distribution').append('<li>' + kx.key + ' <b>(' + nf.format(kx.doc_count) + ')</b></li>');
      }
    }
    $('#submission-count').empty();
    $('#submission-count').append('<div style="font-size:0.8em">Total submissions:</div>' + nf.format(response.hits.total.value));
    $('#export-button').show();
    $('#run-count').empty();
    $('#run-count').append('<div style="font-size:0.8em">Total runs: </div>' + nf.format(response.aggregations.run_count.value));

    const data = response.aggregations.date.buckets;
    let max = 0;
    for (const date of data) {
      if (date.doc_count > max) {
        max = date.doc_count;
      }
    }
    $('#date-chart').empty();
    for (const date of data) {
      const label = date.key_as_string.substring(0, 4);
      if (parseInt(label) > new Date().getFullYear() - 18 && parseInt(label) <= new Date().getFullYear()) {
        const p = (date.doc_count / max) * 100;
        $('#date-chart').append('<li><span style="height:' + p + '%" title="' + label + '"></span></li>');
      }
    }
  }

  async function scroll_page() {
    $('#loading-indicator').show();
    const response = await client.scroll({
      scrollId: window.scroll_id,
      scroll: '10m'
    });
    $('#loading-indicator').hide();
    window.scroll_page += 1;
    update_table(response);
  }

  function update_table(response) {
    $('#page-indicator').empty();
    $('#page-indicator').append('' + ((window.scroll_page - 1) * default_page_size + 1) + ' to ' + ((response.hits.total.value < default_page_size || ((window.scroll_page * default_page_size) > response.hits.total.value)) ? nf.format(response.hits.total.value) : nf.format(window.scroll_page * default_page_size)) + ' of ' + nf.format(response.hits.total.value));
    $('#submissions-table').empty();
    for (const hit of response.hits.hits) {
      let highlights = "";
      if ('highlight' in hit) {
        for (const field of Object.keys(hit.highlight)) {
          for (const match of hit.highlight[field]) {
            highlights += '<div style="margin-bottom: 10px;"><span>' + match + '</span> <span style="opacity: 0.5;font-size:75%">' + field + '</span></div>';
          }
        }
      }
      const st = $('#submissions-table');
      st.append('<tr>');
      st.append('<td style="font-family: monospace;border-bottom: 1px solid #3e4a63;font-size: 14px">' + hit._id + '</td>');
      st.append('<td colspan="4" class="hide-bottomline" style="font-size: 0.8em;font-weight:bold">' + (hit._source.sample1_title ? hit._source.sample1_title : '-') + '; ' + (hit._source.experiment1_title ? hit._source.experiment1_title : '-') + '</td>');
      st.append('</tr>');

      st.append('<tr>');
      st.append('<td></td>');
      st.append('<td>' + highlights + '</td>');
      st.append('<td style="font-family: monospace;font-size: 14px">' + (hit._source.sample_count ? hit._source.sample_count : '-') + '</td>');
      st.append('<td style="font-family: monospace;font-size: 14px">' + (hit._source.run_count ? hit._source.run_count : '-') + '</td>');
      st.append('</tr>');
    }
  }

  async function export_submission_ids() {
    $('#loading-indicator').show();
    const response = await client.search(generate_search($('#query-input').val(), 2000));
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

    let blob = new Blob([all_submission_ids.join('\n') + '\n'], {type: 'text/plain'});
    $('#download-submission-ids').attr('href', (window.URL || window.webkitURL).createObjectURL(blob));
    $('#download-submission-ids')[0].click();
  }

  function update_worldmap(response) {
    const colors = ['#fc8d59', '#ef6548', '#d7301f', '#b30000', '#7f0000', '#9f0000', '#af0000', '#cf0000', '#ef0000', '#ff0000'];
    const canvas = document.getElementById('worldmap');
    const top_crop = 190;
    if (canvas.getContext) {
      const ctx = canvas.getContext('2d');
      const img = new Image();
      img.onload = function () {
        const w = img.width;
        const h = img.height;
        ctx.fillStyle = 'white';
        ctx.fillRect(0, 0, w, h);
        ctx.save();
        ctx.globalAlpha = 0.7;
        ctx.drawImage(img, 0, -top_crop);
        ctx.restore();
        if (response) {
          const  data = response.aggregations.worldmap.buckets;
          if (data.length > 0) {
            const color_step = data[0]['doc_count'] / colors.length;
            for (const area of data) {
              const pos = area['key'].split('/');
              const zoom = parseInt(pos[0]);
              const row_tiles = 2 ** zoom;
              const x = parseInt(pos[1]);
              const y = parseInt(pos[2]);
              ctx.fillStyle = colors[Math.round(area['doc_count'] / color_step) - 1] + 'bb';
              ctx.fillRect((w / row_tiles) * x, (h / row_tiles) * y - top_crop, w / row_tiles, h / row_tiles);
            }
          }
        }
      };
      img.src = 'lib/worldmap.png';
    } else {
    }
  }

  update_worldmap();
});