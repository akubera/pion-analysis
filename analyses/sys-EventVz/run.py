#!/usr/bin/env python3
#
# analysis
#

from pathlib import Path
import yaml

TEMPLATE = """



void
setup_grid(AliAnalysisManager &mgr, TString grid_mode)
{
  auto *alien = new AliAnalysisAlien();
  alien->SetRunMode(grid_mode);
  alien->SetGridOutputDir("output");
  alien->SetGridWorkingDir(workdir);
  alien->SetAliPhysicsVersion("{{ data['aliphysics-version'] }}");
  alien->SetDropToShell(false);
  alien->SetCheckCopy(false);
  alien->SetMaxMergeFiles("{{ data['max-merge-files'] }}");
  alien->SetMaxMergeStages(20);
  alien->SetSplitMaxInputFileNumber("{{ data['max-input-files'] }}");
  alien->SetNrunsPerMaster(30);
  alien->SetMergeViaJDL(true);
  alien->SetTTL({{ data['ttl'] }});
  alien->AddDataFile({{ data['xml-file'] }});

  mgr->SetGridHandler(alien);
}

void
AddTasks()
{
{% for task in analysis['tasks'] %}
    {% if 'args' not in task %} gROOT->Macro("{{ task['macro'] }}");
    {% else %} gROOT->Macro(
R(""{{ task['macro'] }}(
{{ task['args']}}) )""); {% endif %}
{% endfor %}
}

"""

def main():

    __dir__ = Path(__file__).parent

    data = yaml.safe_load((__dir__ / 'data.yaml').read_text())
    analysis = yaml.safe_load((__dir__ / 'analysis.yaml').read_text())
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

    print(local_workdir)


    # print("foo.C(%s)" % analysis['tasks'][-1]['args'].strip())



if __name__ == "__main__":
    exit(main())
