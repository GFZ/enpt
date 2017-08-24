import argparse


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", help="Name")
    parser.add_argument("--parameter_1", help="Input")
    return parser

 
def test(name):
    print(name)


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()    
    test(name=args.name)
