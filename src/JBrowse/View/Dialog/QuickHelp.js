define( [
            'dojo/_base/declare',
            'JBrowse/View/InfoDialog'
        ],
        function(
            declare,
            InfoDialog
        ) {
return declare( InfoDialog, {

    title: "JBrowse Help",

    constructor: function(args) {
        this.browser = args.browser;
    },

    show: function() {
        // make a div containing our help text
        var browserRoot = this.browser.config.browserRoot || "";
        var html = ''
            + '<div class="help_dialog">'
            + '<div class="main" style="float: left; width: 49%;">'

            + '<dl>'
            + '<dt>Moving</dt>'
            + '<dd><ul>'
            + '    <li>Move the view by clicking and dragging in the track area, or by clicking <img height="20px" src="'+browserRoot+'img/slide-left.png"> or <img height="20px"  src="'+browserRoot+'img/slide-right.png"> in the navigation bar, or by pressing the left and right arrow keys.</li>'
            + '    <li>Center the view at a point by clicking on either the track scale bar or overview bar, or by shift-clicking in the track area.</li>'
            + '</ul></dd>'
            + '<dt>Zooming</dt>'
            + '<dd><ul>'
            + '    <li>Zoom in and out by clicking <img height="20px" src="'+browserRoot+'img/zoom-in-1.png"> or <img height="20px"  src="'+browserRoot+'img/zoom-out-1.png"> in the navigation bar, or by pressing the up and down arrow keys while holding down "shift".</li>'
            + '    <li>Select a region and zoom to it ("rubber-band" zoom) by clicking and dragging in the overview or track scale bar, or shift-clicking and dragging in the track area.</li>'
            + '    </ul>'
            + '</dd>'
            + '<dt>Showing Tracks</dt>'
            + '<dd><ul><li>Turn a track on by dragging its track label from the "Available Tracks" area into the genome area, or double-clicking it.</li>'
            + '        <li>Turn a track off by dragging its track label from the genome area back into the "Available Tracks" area.</li>'
            + '    </ul>'
            + '</dd>'
            + '</dl>'
            + '</div>'

            + '<div class="main" style="float: right; width: 49%;">'
            + '<dl>'
            + '<dt>Searching</dt>'
            + '<dd><ul>'
            + '    <li>Jump to a feature or reference sequence by typing its name in the location box and pressing Enter.</li>'
            + '    <li>Jump to a specific region by typing the region into the location box as: <span class="example">ref:start..end</span>.</li>'
            + '    </ul>'
            + '</dd>'
            + '<dt>Example Searches</dt>'
            + '<dd>'
            + '    <dl class="searchexample">'
            + '        <dt>uc0031k.2</dt><dd>searches for the feature named <span class="example">uc0031k.2</span>.</dd>'
            + '        <dt>chr4</dt><dd>jumps to chromosome 4</dd>'
            + '        <dt>chr4:79,500,000..80,000,000</dt><dd>jumps the region on chromosome 4 between 79.5Mb and 80Mb.</dd>'
            + '        <dt>5678</dt><dd>centers the display at base 5,678 on the current sequence</dd>'
            + '    </dl>'
            + '</dd>'
            + '<dt>JBrowse Documentation</dt>'
            + '<dd><ul><li><a target="_blank" href="docs/tutorial/">Quick-start tutorial</a></li>'
            + '        <li><a target="_blank" href="http://gmod.org/wiki/JBrowse">JBrowse Configuration Guide</a></li>'
            + '        <li><a target="_blank" href="docs/config.html"><tt>biodb-to-json.pl</tt> reference</a></li>'
            + '        <li><a target="_blank" href="docs/featureglyphs.html">HTMLFeatures class reference</a></li>'
            + '    </ul>'
            + '</dd>'
            + '</dl>'
            + '</div>'
            + '</div>'
        ;

        this.set( 'content', html );

        this.inherited( arguments );
    }
});
});