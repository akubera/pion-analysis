


class PionAnalysis:

    @dataclass
    class Analysis:
        cuts: dict
        params: dict
        name_prefix: str = 'PionAnalysis_'
        name_suffix: str = ''

        @property
        def name(self):
            charge = self.cuts['track']['charge']
            if charge == 1:
                pair = 'pip'
            elif charge == -1:
                pair = 'pim'
            else:
                raise ValueError(f"Invalid charge {charge!r}")

            cent = '%g_%g' % tuple(self.cuts['event']['cent_range'])

            cfg = self.cfg_obj()

            from ROOT import AliFemtoAnalysisPionPion
            from ROOT import AliFemtoConfigObject
            analysis = AliFemtoAnalysisPionPion(cfg)
            analysis.SetTrackFilter(128)
            outlist = analysis.GetOutputList()
            olist = outlist.At(0)

            cfg = olist.FindObject("AliFemtoConfigObject")
            # ccent = AliFemtoConfigObject.RangeValue_t()
            from hashlib import md5

            m = md5()
            m.update(str(cfg.Stringify()).encode())
            print(m.hexdigest())
            ccent = cfg.pop('event_cut.cent_range')

            # cchrg = AliFemtoConfigObject.IntValue_t()
            cchrg = cfg.pop('track_cut.charge')
            m = md5()
            m.update(str(cfg.Stringify()).encode())

            print(' ->', m.hexdigest())
            # print(cfg.Stringify())
            unique_id = f'{cfg.Hash():016X}_{cent}_{pair}'
            print(cchrg, ccent, unique_id)
            #input()

            return f'{self.name_prefix}{unique_id}{self.name_suffix}'

        def cfg_obj(self):
            from ROOT import AliFemtoConfigObject

            cfg = AliFemtoConfigObject.MapValue_t()
            cfg['event_cut'] = dict_to_config_object(self.cuts['event'])
            cfg['track_cut'] = dict_to_config_object(self.cuts['track'])
            cfg['pair_cut'] = dict_to_config_object(self.cuts['pair'])

            return AliFemtoConfigObject(cfg)

    def __init__(self, cfg):
        self.cfg = dict(cfg)
        self.params = dict(cfg)
        self.multi_cuts = list(self.expand_cuts(self.params.pop('cuts')))
        self.corr_funcs = self.params.pop('corr-fctns')
        # self.params_cfgobj = dict_to_config_object(self.params)

    @staticmethod
    def add_cut_defaults(cuts):

        track = cuts['track']
        if 'charge' not in cuts['track']:
            if 'matrix' not in track:
                track['matrix'] = {"charge": [1, -1]}
            elif 'charge' not in track['matrix']:
                track['matrix']["charge"] = [1, -1]

        event = cuts['event']
        default_centralities = [[0, 10], [20, 30], [40, 50]]
        if 'cent_range' not in event:
            if 'matrix' not in event:
                event['matrix'] = {"cent_range": default_centralities}
            elif 'cent_range' not in event['matrix']:
                event['matrix']["cent_range"] = default_centralities

    @classmethod
    def expand_cuts(cls, cuts):
        """
        Expand multicuts into individual datasets
        """

        cls.add_cut_defaults(cuts)

        event = cuts['event']
        ev_matrix = event.pop('matrix', {})
        events = [{**event, **ev} for ev in expand_matrix(ev_matrix)]

        track = cuts['track']
        tr_matrix = track.pop('matrix', {})
        tracks = [{**track, **tr} for tr in expand_matrix(tr_matrix)]

        pair = cuts['pair']
        pair_matrix = pair.pop('matrix', {})
        pairs = [{**pair, **pr} for pr in expand_matrix(pair_matrix)]

        for ev, tr, pr in _product(events, tracks, pairs):
            yield {"event": ev, "track": tr, "pair": pr}


    @staticmethod
    def read_filter_mask(obj):
        if 'filter-mask' in obj:
            filter_mask = obj['filter-mask']
        elif 'filter-bit' in obj:
            filter_mask = 'BIT(%d)' % int(obj['filter-bit'])
        elif 'filter-bits' in obj:
            filter_mask = '|'.join('BIT(%d)' % b for b in map(int, obj['filter-bit']))
        else:
            raise ValueError("event_reader missing one of these required keys: filter-{bit,bits,mask}")

        return filter_mask

    # def build_event_reader_macro_code(self):
    def evreader_code(self) -> str:

        rdr = self.cfg['event_reader']
        filter_mask = self.read_filter_mask(rdr)

        return '\n'.join([
            f'auto *rdr = new {rdr["_class"]}();',
            f'rdr->SetFilterMask({filter_mask});',
            'rdr->SetReadV0(false);',
            'rdr->SetReadMC(false);',
            'rdr->SetEPVZERO(true);',
            'rdr->SetPrimaryVertexCorrectionTPCPoints(true);',
            'rdr->SetDCAglobalTrack(true);',
            'mgr->SetEventReader(rdr);',
        ])

    def analysis_code(self, analysis) -> str:
        # ev = analysis.cuts['event']
        name = analysis.name

        c = analysis.cfg_obj()

        filter_mask = self.read_filter_mask(self.cfg['event_reader'])

        # c = AliFemtoConfigObject.MapValue_t()
        c.insert("name", name)
        keys = [
            "verbose",
            'events-to-mix',
            'is-simulation',
            "pair-monitors",
            "collection-size-min",
        ]
        for key in keys:
            if key in analysis.params:
                c.insert(key.replace("-", '_'), analysis.params[key])

        configuration = c.Stringify(True)
        return '\n'.join([
            f'auto config = AliFemtoConfigObject::Parse(R"( {configuration} )");',
            f'auto *analysis = new AliFemtoAnalysisPionPion(config);',
            f'analysis->FirstParticleCut()->SetMass(0.13957);',
            f'analysis->AddStanardCutMonitors();',
            f'analysis->SetTrackFilter({filter_mask});',
            'mgr->AddAnalysis(analysis);',
            'SetupCorrelationFunctions(*analysis, 0);',
        ])

    def corrfctn_code(self, cfs) -> str:
        # def
        output = []
        for cf in cfs:
            _cls = cf['class']
            args = stringify_arguments(*cf.get('args', ()))
            kts = cf.get('kt')

            output += ['{', f"auto *cf = new {_cls}({args});"]
            if not kts:
                output += ['analysis.AddCorrFctn(cf);', '}']
            else:

                try:
                    name = kts['name']
                except KeyError:
                    print("Error: correlation function missing name")

                bins = kts['bins']

                output.append(f'auto *kt_binned_cfs = new AliFemtoKtBinnedCorrFunc("{name}", cf);')
                output.extend(f'kt_binned_cfs->AddKtRange({l}, {h}); ' for l, h in bins)
                output.extend([f'analysis.AddCorrFctn(kt_binned_cfs);', '}'])

        return '\n'.join(output + [])

    def get_analyses(self):
        for cut in self.multi_cuts:
            yield self.Analysis(cut, self.params)

    def get_corrfcnts(self):
        yield 0, self.corr_funcs

    def get_template_data(self, filename):
        return {'filename': filename,
                'build_event_reader_macro_code': self.evreader_code,
                'build_analysis': self.analysis_code,
                'build_cf': self.corrfctn_code,
                'get_analyses': self.get_analyses,
                'get_correlation_function_groups': self.get_corrfcnts,
                }

    def dump_config_macro(self, f):
        """
        """
        from pathlib import Path
        from datetime import datetime

        data = self.get_template_data('ConfigFemtoTask.C')

        # for x in self.get_analyses():
        #     # print("X:", x.cuts)
        #     cfg = dict_to_config_object(x.cuts)

        # # for x in self.get_analyses():
        # #     print("Y:", x.name)
        # exit(0)

        if isinstance(f, str):
            f = Path(f)

        if isinstance(f, Path):
            f = f.open('w')

        from jinja2 import Environment, FileSystemLoader

        env = Environment(loader=FileSystemLoader('tmpl'))
        env.globals.update(datetime=datetime)
        env = env.get_template('ConfigFemtoTask.C.jinja')

        txt = env.render(**data)
        f.write(txt)


def expand_matrix(matrix):
    for values in _product(*matrix.values()):
        yield {k: v for k, v in zip(matrix.keys(), values)}


def dict_to_config_object(obj):
    from ROOT import AliFemtoConfigObject
    # builder = AliFemtoConfigObject.BuildMap()
    # builder('key', 'Value')
    cfg = stringify_config_object(obj)
    return AliFemtoConfigObject.Parse(cfg)


def stringify_config_object(cfg):
    def _normalize_keyvalues(kvpair: Tuple[str, Any]):
        key, val = kvpair
        key = key.replace("-", "_")

        if isinstance(val, (tuple, list)) and len(val) >= 2:
            # detect ranges
            if all(isinstance(v, (float, int)) for v in val):
                val = ':'.join('%g' % v for v in val)
        elif isinstance(val, str):
            val = "'%s'" % val
        elif isinstance(val, dict):
            val = _into_configobject_str(val, False)

        return key, val

    def _into_configobject_str(obj: dict, check_class=True) -> str:
        ev_info = dict(obj)
        if check_class:
            import ROOT
            ev_class = ev_info.get("_class")
            if not isinstance(ev_class, str):
                raise TypeError(f"Unexpected type {type(ev_class)}")
            try:
                getattr(ROOT, ev_class)
            except AttributeError:
                print(f"Warning: Could not load class {ev_class} from ROOT")

        items = map(_normalize_keyvalues, ev_info.items())
        strs = (f"{key}: {val}" for key, val in items)
        return '{%s}' % ', '.join(strs)

    return _into_configobject_str(cfg, False)
