import click


class SpecialHelpOrder(click.Group):
    # from https://stackoverflow.com/questions/47972638/how-can-i-define-the-order-of-click-sub-commands-in-help
    def __init__(self, *args, **kwargs):
        self.help_priorities = {}
        super().__init__(*args, **kwargs)

    def get_help(self, ctx):
        self.list_commands = self.list_commands_for_help
        return super().get_help(ctx)

    def list_commands_for_help(self, ctx):
        """reorder the list of commands when listing the help"""
        commands = super().list_commands(ctx)
        return (
            c[1]
            for c in sorted(
                (self.help_priorities.get(command, 1), command) for command in commands
            )
        )

    def command(self, *args, **kwargs):
        """Behaves the same as `click.Group.command()` except capture
        a priority for listing command names in help.
        """
        help_priority = kwargs.pop("help_priority", 1)
        help_priorities = self.help_priorities

        def decorator(f):
            cmd = super(SpecialHelpOrder, self).command(*args, **kwargs)(f)
            help_priorities[cmd.name] = help_priority
            return cmd

        return decorator
