
from dataclasses import dataclass as _dataclass, asdict as _asdict
from pathlib import Path as _Path
from hashlib import md5 as _md5


@_dataclass
class Production:
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

    def xmlfilenames_for_runs(self, runs, parent=None):
        return [self.xmlfilename_for_run(run, parent) for run in runs]

    def xmlfilename_for_run(self, run, parent=None):
        """
        Return unique filename for requested run - no check for valid
        path is performed

        If Parent is provided, result is path relative to parent
        """

        md5 = _md5()
        key = f"{self.data_dir}/{self.prefix}{run}/{self.pattern}"
        md5.update(key.encode())

        xmlfilename = '%d-%s.xml' % (run, md5.hexdigest())

        if parent is not None:
            xmlfilename = _Path(parent) / xmlfilename
        return _Path(xmlfilename)

    def fetch_xml_of_run(self, run, dest=None, sync=True):
        """
        """

        query = self.query_for_run(run)

        def _process_xml(cmd, success, exit_code):
            xml_data = cmd.stdout.strip()
            return xml_data

        import sh

        is_async = not sync

        collection_name = 'collection-%s' % run

        find_proc = sh.alien_find("-x", collection_name, *query,
                                  _done=_process_xml,
                                  _bg=is_async)
        if is_async:
            return find_proc

        return find_proc.wait()

    def xml_of_run(self, run, xmldir=None, sync=True, zip=False):
        """
        Return xml text of query for datafiles in run
        """
        import gzip

        if xmldir is None:
            xmldir = 'xmls'
        xmldir = _Path(xmldir)
        xmldir.mkdir(exist_ok=True, parents=True)

        # check cache locations
        xml_path = self.xmlfilename_for_run(run, xmldir)
        if xml_path.exists():
            return xml_path.read_text()

        zxml_path = xml_path.with_suffix(".xml.gz")
        if zxml_path.exists():
            return gzip.decompress(zxml_path.read_bytes())

        # not in cache - fetch from server
        asyncjob = self.fetch_xml_of_run(run, xmldir, True)
        cmd = asyncjob.wait()
        xml_data = cmd.stdout.strip()

        if zip:
             with gzip.open(zxml_path, "wb") as f:
                 f.write(xml_data)
        else:
            xml_path.write_bytes(xml_data)

        return xml_data.decode()

    def __iter__(self):
        for k, v in _asdict(self).items():
            if isinstance(v, _Path):
                v = str(v)
            yield k, v
