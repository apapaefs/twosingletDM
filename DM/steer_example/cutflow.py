"""Cumulative constraint selections for generic DM scans."""


def stages(rows):
    definitions = [('All evaluated', None), ('Relic density', 'relic'),
                   ('+ direct detection', 'direct_detection'), ('+ gamma lines', 'indirect_detection')]
    if any(row.get('dm_cmb_enabled') for row in rows):
        definitions.append(('+ Planck CMB', 'cmb'))
    verdicts = [True] * len(rows)
    for label, constraint in definitions:
        if constraint:
            for index, row in enumerate(rows):
                excluded = row.get('dm_' + constraint + '_excluded')
                if constraint == 'cmb' and not row.get('dm_cmb_enabled'):
                    excluded = False
                if verdicts[index] is False:
                    continue
                if row.get('dm_calculation_status') != 'success':
                    verdicts[index] = None
                elif excluded is True:
                    verdicts[index] = False
                elif verdicts[index] is None or excluded is None:
                    verdicts[index] = None
        yield {'stage': label, 'passed': sum(v is True for v in verdicts),
               'excluded': sum(v is False for v in verdicts),
               'unassessed': sum(v is None for v in verdicts),
               'indices': [row['index'] for row, verdict in zip(rows, verdicts) if verdict is True]}
