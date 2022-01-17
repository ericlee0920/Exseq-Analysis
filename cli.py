import spatialDE_app
import click


@click.command()
@click.option("-c", "--counts_csv", required=True, type=str, help="csv file for the counts.")
@click.option("-m", "--metadata_csv", required=True, type=str, help="csv file for the metadata with columns center_x, center_y.")
@click.option('--batch-correct/--no-batch-correct', default=False)
def run(**kwargs):
    """Simple program that greets NAME for a total of COUNT times."""
    spatialDE_app.run(**kwargs)

@click.group()
def main():
    pass

main.add_command(run)

if __name__ == '__main__':
    main()
