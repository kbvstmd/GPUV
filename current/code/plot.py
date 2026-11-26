"""Redraw GPUV panels from the plotting inputs released with this snapshot."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
BLUE = '#3B7EA8'
GOLD = '#D4A24C'
TEAL = '#42999A'
GRAY = '#CDD5DA'
REGION_COLORS = ['#AA9271', '#6E9D85', '#4C91A0', '#D27E67', '#617EA5', '#9477A3']


def read(directory, name):
    return pd.read_csv(DATA / directory / f'{name}.tsv', sep='\t')


def save(fig, output, panel):
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f'Fig_{panel}.pdf')
    plt.close(fig)


def diversity(out, figure):
    folder = f'diversity/{figure}'
    points = read(folder, 'coordinates')
    ranks = read(folder, 'rank_summary')
    test = read(folder, 'tests').iloc[0]
    if len(points) != test['n']:
        raise ValueError(f'Figure {figure}: sample count and coordinate count differ')
    title = f"ANOSIM R = {test['ANOSIM_R']:.4f}  •  PERMANOVA R² = {test['PERMANOVA_R2']:.4f}"
    fig, ax = plt.subplots(figsize=(6.1, 4.7), layout='constrained')
    palette = plt.get_cmap('tab20').colors
    for i, (group, values) in enumerate(points.groupby('group', sort=True)):
        color = palette[i]
        ax.scatter(values.PCoA1, values.PCoA2, s=3, alpha=.48, linewidths=0, color=color, label=group)
        center = values[['PCoA1', 'PCoA2']].mean()
        ax.scatter(center.PCoA1, center.PCoA2, s=28, marker='x', linewidths=.9, color=color)
    ax.set(xlabel='PCoA 1', ylabel='PCoA 2', title=title)
    ax.legend(markerscale=2.1, fontsize=6.5, frameon=False, ncols=2)
    save(fig, out, f'{figure}b')

    rank_panel = '2c' if figure == '2' else 'S3b_ranks'
    boxes = [{'label': row.group, 'med': row.median, 'q1': row.q1, 'q3': row.q3,
              'whislo': row.min, 'whishi': row.max} for row in ranks.itertuples()]
    fig, ax = plt.subplots(figsize=(5.6, 4.2), layout='constrained')
    ax.bxp(boxes, orientation='horizontal', showfliers=False, patch_artist=True,
           boxprops={'facecolor': '#D6E6ED'}, medianprops={'color': '#263748'})
    ax.invert_yaxis()
    ax.set_xlabel('Distance rank / total sample pairs')
    save(fig, out, rank_panel)


def candidate_resistance(out):
    folder = 'Fig_5'
    counts = read(folder, 'Panel_A_ARG_class_summary').sort_values('n_references', ascending=True)
    fig, ax = plt.subplots(figsize=(5.8, 4.0), layout='constrained')
    bars = ax.barh(counts.ARG_class, counts.n_references, color=BLUE)
    ax.bar_label(bars, [f'{r.n_references} ({r.percent_of_446:.1f}%)' for r in counts.itertuples()], padding=3, fontsize=7)
    ax.set(xlabel='Reference sequences', xlim=(0, counts.n_references.max() * 1.36))
    save(fig, out, '5a')

    host = read(folder, 'Fig_5d_host_assignment')
    fraction = host.assigned / host.total
    fig, ax = plt.subplots(figsize=(4.8, 2.2), layout='constrained')
    ax.barh(host.candidate_group, fraction, color=BLUE, label='Assigned')
    ax.barh(host.candidate_group, 1 - fraction, left=fraction, color=GRAY, label='Unassigned')
    ax.set(xlabel='Fraction of contigs', xlim=(0, 1))
    ax.legend(frameon=False, ncols=2, loc='lower center', bbox_to_anchor=(.5, 1.02))
    for i, value in enumerate(fraction):
        ax.text(value / 2, i, f'{value:.1%}', ha='center', va='center', color='white')
    save(fig, out, '5d_inset')

    env = read(folder, 'Figure4_environment').sort_values('prevalence', ascending=False)
    x = np.arange(len(env))
    fig, ax = plt.subplots(figsize=(8.4, 3.2), layout='constrained')
    ax.vlines(x - .09, env.wilson_low * 100, env.wilson_high * 100, color=GRAY)
    ax.scatter(x - .09, env.prevalence * 100, color=BLUE, label='Run-weighted')
    ax.scatter(x + .09, env.project_equal_prevalence * 100, color=GOLD, marker='D', label='Project-equal')
    ax.set_xticks(x, [f'{r.name}\n({r.n:,})' for r in env.itertuples()], rotation=35, ha='right')
    ax.set(ylabel='Detection (%)', ylim=(-3, 103))
    ax.legend(frameon=False)
    save(fig, out, '5e')

    locations = read(folder, 'Panel_F_displayed_locations')
    flags = read(folder, 'Location_geography_audit').set_index('city_key').special_geography_flag
    matrix = read(folder, 'Panel_F_heatmap_percent').set_index('ARG_class')[locations.city_key]
    fig, ax = plt.subplots(figsize=(12, 5.2), layout='constrained')
    cmap = LinearSegmentedColormap.from_list('detection', ['#FCFCFA','#D7E8E8','#86BBBD','#337F91','#164965'])
    heatmap = ax.pcolormesh(matrix.to_numpy(), cmap=cmap, vmin=0, vmax=100,
                            edgecolors='white', linewidth=.14)
    ax.invert_yaxis()
    ax.set_yticks(np.arange(len(matrix)) + .5, matrix.index, fontsize=7)
    labels = [f"{r.city_key.split(' | ', 1)[-1]}{'*' if flags.loc[r.city_key] else ''} [{r.country_code}]"
              for r in locations.itertuples()]
    ax.set_xticks(np.arange(len(labels)) + .5, labels, rotation=90, fontsize=5.8)
    for _, block in locations.groupby('region', sort=False):
        start, stop = block.index.min(), block.index.max() + 1
        ax.axvline(start, color='#B9BFC6', linewidth=.6)
        ax.text((start + stop) / 2, -.28, block.region.iloc[0], ha='center', fontsize=6)
    china = locations.index[locations.country_code.eq('CN')]
    y = matrix.index.get_loc('phenicol')
    ax.add_patch(Rectangle((china.min(), y), len(china), 1, fill=False, ec='#C8782B', lw=1.1))
    scale = fig.colorbar(heatmap, ax=ax, orientation='horizontal', shrink=.26, pad=.12)
    scale.solids.set_rasterized(False)
    scale.set_label('Detection (%)')
    save(fig, out, '5f')


def eukaryotic_summary(out):
    folder = 'Fig_6'
    counts = read(folder, 'A_family_counts').sort_values('QC_kept_n', ascending=True)
    fig, ax = plt.subplots(figsize=(5.4, 2.8), layout='constrained')
    ax.barh(counts.family, counts.human_reference_n, label='Human-reference-associated', color='#A60000')
    ax.barh(counts.family, counts.QC_kept_n - counts.human_reference_n, left=counts.human_reference_n,
            label='Other catalogue assignments', color='#DCA3A3')
    ax.set_xlabel('Number of vOTUs')
    ax.legend(frameon=False, fontsize=7)
    save(fig, out, '6a_counts')

    global_counts = read(folder, 'B_global_detection_project_bootstrap').sort_values('prevalence')
    fig, ax = plt.subplots(figsize=(5.4, 3.4), layout='constrained')
    y = np.arange(len(global_counts))
    ax.hlines(y, global_counts.CI_low * 100, global_counts.CI_high * 100, color=GRAY, linewidth=1.6)
    ax.scatter(global_counts.prevalence * 100, y, color=BLUE)
    ax.set_yticks(y, global_counts.endpoint.str.replace('_', ' '))
    ax.set_xlabel('Detection frequency (%)')
    save(fig, out, '6b')

    environment = read(folder, 'C_environment_two_part_summary').query("endpoint == 'Human_reference'")
    environment = environment.sort_values('prevalence')
    fig, axes = plt.subplots(ncols=2, figsize=(7.1, 4.0), sharey=True, layout='constrained')
    ys = np.arange(len(environment))
    ix = environment.n_positive.gt(0)
    axes[0].hlines(ys, environment.CI_low * 100, environment.CI_high * 100, color=GRAY)
    axes[0].scatter(environment.prevalence * 100, ys, color=BLUE)
    axes[0].set_yticks(ys, environment.environment)
    axes[0].set_xlabel('Detection (%)')
    axes[1].hlines(ys[ix], environment.q25_positive[ix] * 100, environment.q75_positive[ix] * 100, color=GRAY)
    axes[1].scatter(environment.median_positive[ix] * 100, ys[ix], color=TEAL)
    axes[1].set(xscale='log', xlabel='Median positive abundance (%)')
    save(fig, out, '6c_summary')

    strata = read(folder, 'D_CMH_full_statistics')
    matrix = strata.pivot(index='endpoint', columns='environment', values='prevalence') * 100
    labels = strata.pivot(index='endpoint', columns='environment', values='positive_only_annotation')
    fig, ax = plt.subplots(figsize=(9.6, 4.0), layout='constrained')
    mesh = ax.pcolormesh(matrix.to_numpy(), vmin=0, vmax=50, cmap='Reds')
    for i, endpoint in enumerate(matrix.index):
        for j, environment_name in enumerate(matrix.columns):
            val = matrix.loc[endpoint, environment_name]
            tag = labels.loc[endpoint, environment_name]
            ax.text(j+.5, i+.5, f'{val:.1f}\n{tag}' if pd.notna(tag) else f'{val:.1f}',
                    va='center', ha='center', fontsize=6.2, color='white' if val > 32 else '#252525')
    ax.set_xticks(np.arange(len(matrix.columns)) + .5, matrix.columns, rotation=45, ha='right')
    ax.set_yticks(np.arange(len(matrix.index)) + .5, matrix.index.str.replace('_', ' '))
    ax.invert_yaxis()
    scale = fig.colorbar(mesh, ax=ax, label='Detection (%)')
    scale.solids.set_rasterized(False)
    save(fig, out, '6d')


def host_assignment(out):
    for panel, colors in [('S4a', [GRAY, BLUE, GOLD]), ('S4b', [BLUE, TEAL, '#7466A3', GRAY])]:
        summary = read('Fig_S4', f'Fig_{panel}')
        fig, ax = plt.subplots(figsize=(5.0, 2.7), layout='constrained')
        wedges, _ = ax.pie(summary.value, colors=colors, startangle=90, counterclock=False,
                           wedgeprops={'width': .34, 'edgecolor': 'white'})
        ax.legend(wedges, summary.category + '  ' + summary.label, loc='center left',
                  bbox_to_anchor=(.96, .5), frameon=False)
        save(fig, out, panel)


def carrier_counts(out):
    folder = 'Fig_S7/carrier_counts'
    counts = read(folder, 'count_histograms')
    tests = read(folder, 'overall_tests')
    factors = {
        'environment': ('g', ['WWTP','Indoor','Sludge','Biofilm','Hospital','Transport hubs',
                              'Freshwater','Indoor air','Outdoor','Bioreactor','Air (unspecified)','Coast','Soil']),
        'continent': ('h', ['Africa','Europe','Asia','South America','North America','Oceania']),
    }
    for factor, (letter, groups) in factors.items():
        selected = tests.loc[tests.cohort.eq('all') & tests.factor.eq(factor) & tests.model.eq('unadjusted')].iloc[0]
        histogram = counts.loc[counts.factor.eq(factor)]
        values = []
        for group in groups:
            rows = histogram.loc[histogram.group.eq(group)]
            values.append(np.repeat(rows['count'], rows.frequency))
        if sum(map(len, values)) != selected['n'] or any(len(group) == 0 for group in values):
            raise ValueError(f'{factor}: histogram and plotted cohort differ')
        fig, ax = plt.subplots(figsize=(6.5 if letter == 'g' else 4.0, 3.4), layout='constrained')
        ax.boxplot(values, tick_labels=groups, patch_artist=True,
                   boxprops={'facecolor': '#E2ECEC', 'edgecolor': '#557F83'},
                   medianprops={'color': '#263748'}, flierprops={'markersize': 2, 'markeredgecolor': '#557F83'})
        ax.tick_params(axis='x', rotation=50)
        for label in ax.get_xticklabels():
            label.set_ha('right')
        ax.set(ylabel='Detected carrier vOTUs per run', ylim=(0, None), title=f"Overall P = {selected['p']:.4f}")
        save(fig, out, f'S7{letter}')


def main():
    plots = {'2': lambda out: diversity(out, '2'), 'S3': lambda out: diversity(out, 'S3'),
             '5': candidate_resistance, '6': eukaryotic_summary, 'S4': host_assignment, 'S7': carrier_counts}
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fig', choices=plots, nargs='+', default=list(plots))
    parser.add_argument('--out', type=Path, default=ROOT / 'pdf' / 'redrawn')
    args = parser.parse_args()
    font_manager.findfont('Arial', fallback_to_default=False)
    plt.rcParams.update({'font.family':'Arial', 'font.size':9, 'pdf.fonttype':42,
                         'pdf.use14corefonts':False, 'axes.spines.top':False, 'axes.spines.right':False})
    for key in args.fig:
        plots[key](args.out)


if __name__ == '__main__':
    main()
