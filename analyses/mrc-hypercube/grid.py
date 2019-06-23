#!/usr/bin/env python
#
# grid.py
#


import sys

def argparser():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("task",
                        choices=('submit','merge'))
    parser.add_argument("job", default=None)
    return parser


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = argparser().parse_args(argv)
    if args.task == 'submit':
        return submit(args)


def submit(args):
    print("submitting job:", args.job)


if __name__ == "__main__":
    exit(main())
