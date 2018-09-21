def write_file(filename, contents):
    outfile = open(filename, 'w')
    outfile.write(contents)
    outfile.close()


def get_opts(argv):
    opts = {}
    while argv:
        if argv[0][0] == '-':
            opts[argv[0]] = argv[1]
        argv = argv[1:]
    return opts