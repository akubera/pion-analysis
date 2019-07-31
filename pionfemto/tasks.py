#
# pionfemto/tasks.py
#

# from femtoanalysis.task import Task


class Task:
    """
    """

    @classmethod
    def From(cls, obj) -> 'Task':
        keyset = set(obj.keys())
        if keyset == {'macro'}:
            return MacroTask(obj)
        elif 'analysis' in keyset:
            return AnalysisTask(obj)


class MacroTask(Task):
    """ Task composed of a single call to a macro """

    def __init__(self, obj):
        self.macro = obj['macro']


class AnalysisTask(Task):

    def __init__(self, obj):

        from . import analyses

        analysis_classname = obj['analysis']

        analysis_cls = getattr(analyses, analysis_classname, None)
        if analysis_classname:
            analysis = analysis_cls(obj)

        # pprint(femtoanalysia.analysis[])
        return analysis
