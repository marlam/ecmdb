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
#include <cstdlib>

#include "msg.h"
#include "opt.h"
#include "str.h"

#include "cmds.h"


// The code to manage and start the commands was taken from cvtool-0.2.6.

/*
 * The command functions. All live in their own .cpp file, except for the
 * trivial help and version commands. Some are builtin, and some are loaded
 * as modules (notably those that depend on external libraries).
 * Add the aproppriate lines for each command in this section.
 * The commands must appear in ascending order according to strcmp(), because we
 * do a binary search on the command name.
 */

typedef struct
{
    const char *name;
    int (*cmd)(int argc, char *argv[]);
    void (*cmd_print_help)(void);
} cmd_t;

#define CMD_DECL(FNBASE) \
    extern "C" int ecmdb_ ## FNBASE (int argc, char *argv[]); \
    extern "C" void ecmdb_ ## FNBASE ## _help (void);

#define CMD(NAME, FNBASE) { NAME, ecmdb_ ## FNBASE, ecmdb_ ## FNBASE ## _help }

CMD_DECL(add)
CMD_DECL(commit)
CMD_DECL(create)
CMD_DECL(filename)
CMD_DECL(finalize)
CMD_DECL(help)
CMD_DECL(resolutions)
CMD_DECL(version)

static cmd_t cmds[] =
{
    CMD("add",         add        ),
    CMD("commit",      commit     ),
    CMD("create",      create     ),
    CMD("filename",    filename   ),
    CMD("finalize",    finalize   ),
    CMD("help",        help       ),
    CMD("resolutions", resolutions),
    CMD("version",     version    ),
};


/*
 * If you just want to add a command, there's no need to change anything below
 * this line.
 */

int cmd_count()
{
    return (sizeof(cmds) / sizeof(cmds[0]));
}

const char *cmd_name(int cmd_index)
{
    return cmds[cmd_index].name;
}

static int cmd_strcmp(const void *a, const void *b)
{
    const cmd_t *c1 = static_cast<const cmd_t *>(a);
    const cmd_t *c2 = static_cast<const cmd_t *>(b);
    return strcmp(c1->name, c2->name);
}

int cmd_find(const char *cmd)
{
    cmd_t *p;
    cmd_t key = { cmd, NULL, NULL };

    p = static_cast<cmd_t *>(bsearch(
                static_cast<void *>(&key),
                static_cast<void *>(cmds),
                cmd_count(),
                sizeof(cmd_t),
                cmd_strcmp));
    int cmd_index = (p ? p - cmds : -1);
    return cmd_index;
}

void cmd_run_help(int cmd_index)
{
    cmds[cmd_index].cmd_print_help();
}

int cmd_run(int cmd_index, int argc, char *argv[])
{
    return cmds[cmd_index].cmd(argc, argv);
}
