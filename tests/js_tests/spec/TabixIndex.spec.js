require({
            packages: [{ name: 'jDataView', main: 'jdataview', location: '../jDataView/src' }]
        },
        [
            'JBrowse/Model/TabixIndex',
            'JBrowse/Model/XHRBlob',
            'JBrowse/Model/BGZip/BGZBlob'
        ],function( TabixIndex, XHRBlob, BGZBlob ) {

describe( "TabixIndex", function() {

              var ti;
              beforeEach(function() {
                  ti = new TabixIndex( new BGZBlob( new XHRBlob( '../../sample_data/raw/volvox/volvox.test.vcf.gz.tbi' ) ) );
              });

              it('constructs', function() {
                     expect( ti ).toBeTruthy();
                 });
              it( 'loads', function() {
                      var loaded;
                      waitsFor( function() { return loaded; }, 1000 );
                      ti.load().then( function() {
                          loaded = true;
                      });
                      runs( function() {
                          expect( ti.columnNumbers.start ).toEqual( 2 );
                          expect( ti.columnNumbers.ref ).toEqual( 1 );
                          expect( ti.columnNumbers.end ).toEqual( 0 );
                          console.log( ti );
                          var blocks = ti.blocksForRange( 'ctgA', 1, 4000 );
                          expect( blocks.length ).toEqual( 1 );
                          expect( blocks[0].minv.block ).toEqual( 0 );
                          expect( blocks[0].minv.offset ).toEqual( 10501 );
                          console.log( blocks );
                      });
              });
});
});
