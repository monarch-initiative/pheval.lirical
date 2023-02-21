import click

from pheval_lirical.post_process.post_process_results_format import post_process
from pheval_lirical.prepare.prepare_commands import prepare_commands_command


def main_():
    pass


@click.group()
def main():
    """Lirical runner."""
    pass


main.add_command(prepare_commands_command)
main.add_command(post_process)

if __name__ == "__main__":
    main()
