#!/usr/bin/env python3
#
# analysis
#

import sys
from pathlib import Path
import json
import yaml


def arg_parser():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--data',
                        default='data.yaml')
    parser.add_argument('--analysis',
                        default='analysis.yaml')
    return parser


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    args = arg_parser().parse_args(argv)

    __dir__ = Path(__file__).absolute().parent

    data = yaml.safe_load((__dir__ / args.data).read_text())
    analysis = yaml.safe_load((__dir__ / args.analysis).read_text())
    macro_template = __dir__.parent.parent / 'templates' / 'RunAnalysis.C.j2'

    from alimaster.datasets import Production
    from jinja2 import Template

    macro_txt = macro_template.read_text()
    t = Template(macro_txt)

    production = Production.From(data['production'])

    vals = {
        'aliphysics_version': analysis['aliphysics-version'],
        'TTL': analysis.get('ttl', '8 * 60 * 60'),
        'production': production,
        'analysis': analysis,
    }

    from ROOT import TGrid

    grid = TGrid.Connect("alien://")

    from datetime import datetime
    tag = f'{datetime.now():%y%m%d%H%M%S}'

    local_workdir = __dir__ / "gridresults" / f'job-{tag}'
    local_workdir.mkdir(parents=True)

    for i, (dset_name, runs) in enumerate(data['datasets'].items(), 1):
        hashname = production.get_hashname(runs)

        xmlfilename = production.build_xml_dataset(runs, grid=grid)
        workdir = 'job-%s-%02d' % (tag, i)
        result_filename = 'AnalysisResults-%s-%s.root' % (tag, dset_name)
        dest_wrkdir = local_workdir / dset_name
        dest_wrkdir.mkdir()

        txt = t.render(**vals,
                       xmlfile=Path(grid.GetHomeDirectory()) / 'xml' / xmlfilename,
                       workdir=workdir,
                       output_filename=result_filename)

        (dest_wrkdir / 'RunAnalysis.C').write_text(txt)

    (local_workdir / 'settings.json').write_text(json.dumps({'analysis': analysis, 'data': data}, indent=1))

    print(local_workdir)


if __name__ == "__main__":
    exit(main())
