/*
 * Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012
 * Martin Lambers <marlam@marlam.de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#include <cstring>

#if W32
# include <stdlib.h>
# include <io.h>
# include <fcntl.h>
# include <strings.h>
#endif

#include "msg.h"
#include "opt.h"
#include "dbg.h"

#include "cmds.h"


char* program_name = NULL;
FILE* ecmdb_stdin = NULL;
FILE* ecmdb_stdout = NULL;

extern "C" void ecmdb_version_help(void)
{
    msg::req_txt(
            "version\n"
            "\n"
            "Print version information.");
}

extern "C" int ecmdb_version(int argc, char* argv[])
{
    std::vector<opt::option*> options;
    opt::info help("help", '\0', opt::optional);
    options.push_back(&help);
    std::vector<std::string> arguments;
    if (!opt::parse(argc, argv, options, 0, 0, arguments)) {
        return 1;
    }
    if (help.value()) {
        ecmdb_version_help();
        return 0;
    }
    msg::req_txt("%s (from package %s) version %s on %s\n"
            "Copyright (C) 2012  Martin Lambers <marlam@marlam.de>.\n"
            "All rights reserved.\n"
            "There is NO WARRANTY, to the extent permitted by law.",
            program_name, PACKAGE_NAME, VERSION, PLATFORM);
    return 0;
}

extern "C" void ecmdb_help_help(void)
{
    msg::req_txt(
            "help [<command>]\n"
            "\n"
            "Print general or command specific help.");
}

extern "C" int ecmdb_help(int argc, char *argv[])
{
    std::vector<opt::option*> options;
    opt::info help("help", '\0', opt::optional);
    options.push_back(&help);
    std::vector<std::string> arguments;
    if (!opt::parse(argc, argv, options, 0, 1, arguments)) {
        return 1;
    }
    if (help.value()) {
        ecmdb_help_help();
        return 0;
    }
    if (arguments.size() == 0) {
        msg::req_txt(
                "Usage: %s [-q|--quiet] [-v|--verbose] <command> [argument...]\n",
                program_name);
        for (int j = 0; j < cmd_count(); j++) {
            msg::req("%s", cmd_name(j));
        }
        msg::req_txt(
                "\n"
                "Use \"%s help <command>\" for command specific help.\n"
                "Report bugs to <%s>.", program_name, PACKAGE_BUGREPORT);
        return 0;
    } else {
        int cmd_index = cmd_find(argv[1]);
        if (cmd_index < 0) {
            msg::err("command unknown: %s", argv[1]);
            return 1;
        } else {
            cmd_run_help(cmd_index);
            return 0;
        }
    }
}

int main(int argc, char *argv[])
{
    int exitcode = 0;
#if W32
    _setmode(_fileno(stdin), _O_BINARY);
    _setmode(_fileno(stdout), _O_BINARY);
    _fmode = _O_BINARY;
    setbuf(stderr, NULL);
    program_name = strrchr(argv[0], '\\');
    program_name = program_name ? program_name + 1 : argv[0];
    size_t program_name_len = strlen(program_name);
    if (program_name_len > 4 && strcasecmp(program_name + program_name_len - 4, ".exe") == 0)
    {
        program_name = strdup(program_name);
        program_name[program_name_len - 4] = '\0';
    }
#else
    program_name = std::strrchr(argv[0], '/');
    program_name = program_name ? program_name + 1 : argv[0];
#endif
    msg::set_level(msg::INF);
    msg::set_program_name(program_name);
    msg::set_columns_from_env();
    dbg::init_crashhandler();

    if (argc < 2) {
        char help[] = "help";
        char *my_argv[] = { help, NULL };
        ecmdb_help(1, my_argv);
        exitcode = 1;
    } else if (argc == 2 && strcmp(argv[1], "--help") == 0) {
        char help[] = "help";
        char *my_argv[] = { help, NULL };
        exitcode = ecmdb_help(1, my_argv);
    } else if (argc == 2 && strcmp(argv[1], "--version") == 0) {
        char version[] = "version";
        char *my_argv[] = { version, NULL };
        exitcode = ecmdb_version(1, my_argv);
    } else {
        int argv_cmd_index = 1;
        if (argc > argv_cmd_index + 1 && (strcmp(argv[argv_cmd_index], "-q") == 0
                    || strcmp(argv[argv_cmd_index], "--quiet") == 0)) {
            argv_cmd_index++;
            msg::set_level(msg::WRN);
        }
        if (argc > argv_cmd_index + 1 && (strcmp(argv[argv_cmd_index], "-v") == 0
                    || strcmp(argv[argv_cmd_index], "--verbose") == 0)) {
            argv_cmd_index++;
            msg::set_level(msg::DBG);
        }
        int cmd_index = cmd_find(argv[argv_cmd_index]);
        if (cmd_index < 0) {
            msg::err("command unknown: %s", argv[argv_cmd_index]);
            exitcode = 1;
        } else {
            msg::set_program_name(msg::program_name() + " " + argv[argv_cmd_index]);
            ecmdb_stdin = stdin;
            ecmdb_stdout = stdout;
            exitcode = cmd_run(cmd_index, argc - argv_cmd_index, &(argv[argv_cmd_index]));
        }
    }
    return exitcode;
}
