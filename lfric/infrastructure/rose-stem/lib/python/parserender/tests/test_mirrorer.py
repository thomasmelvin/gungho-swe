##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Ensure all functionality relating to mirroring is working.

TODO: These tests are extremely incomplete.
"""
from pathlib import Path

from ..mirrorer import Mirrorer


class TestMirrorer:
    """
    Ensure the mirroring class is functioning correctly.

    TODO: One recent changes are reflected in the testing.
    """
    def test_entities(self, tmp_path: Path):
        """
        Makes sure we correctly handle character entities.
        """
        inject_text = """<html>
  <head>
    <title>Test character entities</title>
  </head>
  <body>
    <ul>
      <li>&nbsp;</li>
      <li>&quot;</li>
      <li>&cent;</li>
      <li>&pound;</li>
      <li>&sect;</li>
      <li>&copy;</li>
      <li>&laquo;</li>
      <li>&raquo;</li>
      <li>&reg;</li>
      <li>&deg;</li>
      <li>&plusmn;</li>
      <li>&para;</li>
      <li>&middot;</li>
      <li>&frac12;</li>
      <li>&ndash;</li>
      <li>&mdash;</li>
      <li>&lsquo;</li>
      <li>&rsquo;</li>
      <li>&sbquo;</li>
      <li>&ldquo;</li>
      <li>&rdquo;</li>
      <li>&bdquo;</li>
      <li>&dagger;</li>
      <li>&Dagger;</li>
      <li>&bull;</li>
      <li>&hellip;</li>
      <li>&prime;</li>
      <li>&Prime;</li>
      <li>&euro;</li>
      <li>&trade;</li>
      <li>&asymp;</li>
      <li>&ne;</li>
      <li>&le;</li>
      <li>&ge;</li>
      <li>&lt;</li>
      <li>&gt;</li>
      <li>&amp;</li>
    </ul>
    <p>Multiple &lsquo;in one line&rsquo;.</p>
  </body>
</html>"""
        expect_text = """<html>
  <head>
    <title>Test character entities</title>
  </head>
  <body>
    <ul>
      <li>&#160;</li>
      <li>"</li>
      <li>&#162;</li>
      <li>&#163;</li>
      <li>&#167;</li>
      <li>&#169;</li>
      <li>&#171;</li>
      <li>&#187;</li>
      <li>&#174;</li>
      <li>&#176;</li>
      <li>&#177;</li>
      <li>&#182;</li>
      <li>&#183;</li>
      <li>&#188;</li>
      <li>&#8211;</li>
      <li>&#8212;</li>
      <li>&#8216;</li>
      <li>&#8217;</li>
      <li>&#8218;</li>
      <li>&#8220;</li>
      <li>&#8221;</li>
      <li>&#8222;</li>
      <li>&#8224;</li>
      <li>&#8225;</li>
      <li>&#8226;</li>
      <li>&#8230;</li>
      <li>&#8242;</li>
      <li>&#8243;</li>
      <li>&#8364;</li>
      <li>&#8482;</li>
      <li>&#8776;</li>
      <li>&#8800;</li>
      <li>&#8804;</li>
      <li>&#8805;</li>
      <li>&lt;</li>
      <li>&gt;</li>
      <li>&amp;</li>
    </ul>
    <p>Multiple &#8216;in one line&#8217;.</p>
  </body>
</html>"""
        test_dirname = tmp_path / 'in'
        test_dirname.mkdir()
        test_filename = test_dirname / 'test.html'
        test_filename.write_text(inject_text)
        generated_dirname = tmp_path / 'out'
        test_unit = Mirrorer(generated_dirname.as_uri())
        test_unit.mirror(test_dirname)
        generated_filename = generated_dirname / 'test.html'
        assert generated_filename.read_text() == expect_text
