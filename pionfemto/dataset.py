
from dataclasses import dataclass as _dataclass, asdict as _asdict
from pathlib import Path as _Path
from hashlib import md5 as _md5
import subprocess as sp


@_dataclass
class Production:
    name: str
    prefix: str
    data_dir: str
    pattern: str
    rootfile: str
    archive: str
    is_mc: bool

    def query_for_run(self, run):
        path = '%s/%s%d' % (self.data_dir, self.prefix, run)
        path += '/' + self.pattern.rstrip("*")
        filename = self.rootfile
        return path, filename

    def xmlfilenames_for_runs(self, runs, parent=None, zip=False):
        yield from (self.xmlfilename_for_run(run, parent) for run in runs)

    def xmlfilename_for_run(self, run, parent=None, zip=False):
        """
        Return unique filename for requested run - no check for valid
        path is performed

        If Parent is provided, result is path relative to parent
        """

        md5 = _md5()
        key = f"{self.data_dir}/{self.prefix}{run}/{self.pattern}"
        md5.update(key.encode())

        xmlfilename = _Path('%d-%s.xml' % (run, md5.hexdigest()))
        if zip:
            xmlfilename = xmlfilename.with_suffix('.xml.gz')

        if parent is not None:
            xmlfilename = _Path(parent) / xmlfilename
        return _Path(xmlfilename)

    def fetch_xml_of_run(self, run, sync=True) -> bytes:
        """
        """
        from sys import stderr
        import subprocess as sp
        from xml.etree import ElementTree

        query = self.query_for_run(run)

        collection_name = 'collection-%s' % run
        find_proc = sp.Popen(['alien_find', '-x', collection_name, *query],
                             stdout=sp.PIPE, stderr=sp.PIPE)

        xml_out, xml_err = find_proc.communicate()

        # "validate" the xml
        try:
            ElementTree.fromstring(xml_out)
        except ElementTree.ParseError:
            print(xml_out.decode().strip(), file=stderr)
            exit(1)

        return xml_out.strip()

    def xml_of_run(self, run, xmldir=None, zip=False, return_path=False) -> str:
        """
        Return xml text of query for datafiles in run
        """
        import gzip

        xmldir = _Path(xmldir or 'xmls')
        xmldir.mkdir(exist_ok=True, parents=True)

        # check cache locations
        xml_path = self.xmlfilename_for_run(run, xmldir)
        if xml_path.exists():
            return xml_path.read_text()

        zxml_path = xml_path.with_suffix(".xml.gz")
        if zxml_path.exists():
            return gzip.decompress(zxml_path.read_bytes()).decode()

        # not in cache - fetch from server
        xml_data = self.fetch_xml_of_run(run)

        if zip:
             with gzip.open(zxml_path, "wb") as f:
                 f.write(xml_data)
        else:
            xml_path.write_bytes(xml_data)

        if return_path:
            return zxml_path if zip else xml_path

        return xml_data.decode()

    def __iter__(self):
        for k, v in _asdict(self).items():
            if isinstance(v, _Path):
                v = str(v)
            yield k, v
