import argparse
import sys


def run_pca(args):
    from src.db.load_data import load_and_prepare_data
    from src.analysis.pca import run_pca, plot_pca

    df = load_and_prepare_data()
    result_df, explained = run_pca(df)
    plot_pca(result_df, explained,
             output_path="data/outputs/pca_plot.png")


def run_outliers(args):
    from src.db.load_data import load_and_prepare_data
    from src.analysis.outlier_detection import (
        detect_outliers, print_outlier_report, plot_outliers
    )
    import yaml

    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    threshold = config['analysis']['outlier_threshold']

    df = load_and_prepare_data()
    outlier_df = detect_outliers(df, threshold=threshold)
    print_outlier_report(outlier_df)
    plot_outliers(df, outlier_df,
                  output_path="data/outputs/outlier_plot.png")


def run_compare(args):
    from src.db.load_data import load_and_prepare_data
    from src.analysis.comparison import (
        group_by_taxonomy, compare_metrics,
        print_comparison_report, plot_comparison
    )

    df = load_and_prepare_data()
    rank = args.rank

    group_by_taxonomy(df, rank=rank)
    summary = compare_metrics(df, rank=rank)
    print_comparison_report(summary, rank)
    plot_comparison(
        df, rank=rank,
        metric=args.metric,
        output_path=f"data/outputs/comparison_{rank}.png"
    )


def run_report(args):
    from src.services.report_service import generate_report
    generate_report(assembly_id=args.assembly_id)


def main():
    parser = argparse.ArgumentParser(
        prog='genome-metrics',
        description='Genome Assembly Metrics Analysis Tool'
    )

    subparsers = parser.add_subparsers(
        title='commands',
        dest='command',
        metavar='<command>'
    )

    # PCA command
    subparsers.add_parser(
        'pca',
        help='Run PCA on genome assembly metrics'
    )

    # Outliers command
    subparsers.add_parser(
        'outliers',
        help='Detect outlier genomes using z-score'
    )

    # Compare command
    compare_parser = subparsers.add_parser(
        'compare',
        help='Compare metrics across taxonomy groups'
    )
    compare_parser.add_argument(
        '--rank',
        default='family',
        choices=['species', 'genus', 'family', 'order_name'],
        help='Taxonomy rank for grouping (default: family)'
    )
    compare_parser.add_argument(
        '--metric',
        default='contig_N50',
        help='Metric to visualize (default: contig_N50)'
    )

    # Report command
    report_parser = subparsers.add_parser(
        'report',
        help='Generate per-genome metrics report'
    )
    report_parser.add_argument(
        '--assembly-id',
        type=int,
        required=True,
        help='Assembly ID to generate report for'
    )

    # Print help if no command given
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Route to correct function
    commands = {
        'pca': run_pca,
        'outliers': run_outliers,
        'compare': run_compare,
        'report': run_report
    }

    commands[args.command](args)


if __name__ == "__main__":
    main()