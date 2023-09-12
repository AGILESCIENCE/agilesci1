////////////////////////////////////////////////////////////////////////////////
// DESCRIPTION
//       AGILE Science Tools
//       AG_indexgen
//       Contributors: 
//       Author: Andrea Zoli (INAF/IASF Bologna)
//
// NOTICE
//       Any information contained in this software
//       is property of the AGILE TEAM and is strictly
//       private and confidential.
//       Copyright (C) 2005-2019 AGILE Team. All rights reserved.
/*
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
////////////////////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <zlib.h>
#include <limits.h>

#define CHUNK_SIZE 64*1024

int main(int argc, char *argv[])
{
    char *log_dir, *out_file;
    char *type;
    DIR *dp;
    struct dirent *ep;
    double tstart, tstop;
    char *out_str;
    FILE *fd, *fdo;
    char buf[PATH_MAX+1], logfile[PATH_MAX+1];
    z_stream strm;
    int ret;
    char in[CHUNK_SIZE], out[CHUNK_SIZE];
    char *ptr;

    if (argc != 4)
    {
        printf("USAGE: %s <log_dir> <type> <out_file>\n", argv[0]);
        return EXIT_FAILURE;
    }
    if (argv[1][0] != '/')
    {
        fprintf(stderr, "%s: error: expected an absolute directory name for --prefix: %s\n", argv[0], argv[1]);
        return EXIT_FAILURE;
    }
    if (argv[3][0] != '/')
    {
        fprintf(stderr, "%s: error: expected an absolute directory name for --prefix: %s\n", argv[0], argv[3]);
        return EXIT_FAILURE;
    }
    log_dir = argv[1];
    type = argv[2];
    out_file = argv[3];
    dp = opendir(log_dir);
    if (!dp)
    {
        fprintf(stderr, "%s: error: error opening directory %s\n", argv[0], log_dir);
        return EXIT_FAILURE;
    }
    out_str = (char*) malloc (52428800);
    if(!out_str)
    {
        fprintf(stderr, "%s: error: allocation error.\n", argv[0]);
        return EXIT_FAILURE;
    }
    out_str[0] = '\0';
    while( (ep = readdir(dp)) )
    {
        sprintf(logfile, "%s/%s", log_dir, ep->d_name);
        if(strstr(logfile, ".gz") == NULL)
            continue;

        fd = fopen(logfile, "rb");
        if (!fd)
        {
            fprintf(stderr, "%s: warn: cannot open file %s.\n", argv[0], logfile);
            continue;
        }

        strm.zalloc = Z_NULL;
        strm.zfree = Z_NULL;
        strm.opaque = Z_NULL;
        strm.avail_in = fread(in, 1, CHUNK_SIZE, fd);
        strm.next_in = (unsigned char*)in;
        strm.avail_out = CHUNK_SIZE;
        strm.next_out = (unsigned char*)out;
        ret = inflateInit2(&strm, 47);
        if (ret != Z_OK)
        {
            fprintf(stderr, "%s: err: cannot initialize inflate\n", argv[0]);
            fclose(fd);
            closedir(dp);
            free(out_str);
            return EXIT_FAILURE;
        }
        ret = inflate(&strm, Z_NO_FLUSH);

        tstart = -1.;
        tstop = -1.;
        ptr = strtok(out, " ");
        while (ptr != NULL && (tstart < 0. || tstop < 0.))
        {
            if(strcmp(ptr, "TSTART") == 0)
            {
                ptr = strtok(NULL, " "); /* = */
                ptr = strtok(NULL, " "); /* value */
                tstart = atof(ptr);
            }
            if(strcmp(ptr, "TSTOP") == 0)
            {
                ptr = strtok(NULL, " "); /* = */
                ptr = strtok(NULL, " "); /* value */
                tstop = atof(ptr);
            }
            ptr = strtok(NULL, " ");
        }

        inflateEnd(&strm);
        fclose(fd);

        sprintf(out_str+strlen(out_str), "%s %lf %lf %s\n", logfile, tstart, tstop, type);
    }

    if(strlen(out_str) > 0)
        out_str[strlen(out_str)-1] = '\0';

    fdo = fopen(out_file, "w");
    if (!fdo)
    {
        fprintf(stderr, "%s: err: cannot open outout file %s.\n", argv[0], out_file);
        closedir(dp);
        free(out_str);
        return EXIT_FAILURE;
    }
    /*printf("%s\n", out_str);*/
    fwrite(out_str, 1, strlen(out_str), fdo);
    fclose(fdo);

    closedir(dp);
    free(out_str);

    return EXIT_SUCCESS;
}
