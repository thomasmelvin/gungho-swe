"""
PSyclone transformation script.
"""


def trans(psy):
    """
    Duplicates the body of the first loop twice leaving three copies.
    """
    # Loop over all of the Invokes in the PSy object
    for invoke in psy.invokes.invoke_list:

        print(f"Transforming invoke '{invoke.name}' ...")
        schedule = invoke.schedule

        loop = schedule.loops()[0]
        loop_schedule = loop.loop_body
        new_node = loop_schedule[0].copy()
        loop_schedule.addchild(new_node)
        add_node = loop_schedule[0].copy()
        loop_schedule.addchild(add_node)

        # Take a look at what we've done
        schedule.view()

    return psy
