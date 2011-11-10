function on_load() {
    'use strict';
    var n_entries = 0;

    function load_list(data) {
        var columns = data.columns.map(function (c) { return c.name; });
        data.columns.forEach(function (c) {
            $('<th/>').html(c.name).css('width', c.width).appendTo('#header');
        });

        data.entries.forEach(function (c) {
            n_entries += 1;
            var row = $('<tr/>').attr('class', 'entry');
            columns.forEach(function (col) {
                $('<td/>').html(c.fields[col]).appendTo(row);
            });
            row.click(function () { $('#content-pane').load(c.content_url); });
            row.dblclick(function () { $('#content-pane').load(c.content_url); hide_panel(); });
            row.appendTo('#list-table');
        });
        show_panel();
    }


    function show_panel() {
        var offset = (n_entries+2) * 1.5 + 0.7;
        $('#list-table').css('display', 'table');
        $('#list-pane').css('height', offset + 'em');
        $('#content-pane').css('top', offset + 'em');
        $('#hide').html('&#x25B2;Hide');
        $('#hide').click(hide_panel);
    }

    function hide_panel() {
        $('#list-table').css('display', 'none');
        $('#list-pane').css('height', '1.5em');
        $('#content-pane').css('top', '1.5em');
        $('#hide').html('&#x25BC;Show');
        $('#hide').click(show_panel);
    }

    $(document).ready(function () {
        $.getJSON('http://localhost/~ross/web_interface/list.json', load_list);
    });
}

on_load();
