// ############################################################################
// STing
// Copyright (C) 2019 Georgia Institute of Technology
// STing is freely available to personal, academic and non-profit use only. You
// cannot redistribute STing to other users. No liability for software usage is
// assumed. For more information on how to obtain a commercial license please 
// conact Lavanya Rishishwar <lavanya.rishishwar@gatech.edu>
// ============================================================================ 
// (c) Georgia Institute of Technology 2019
// Author:  Hector F. Espitia-Navarro
//          hspitia@gatech.edu
//          School of Biological Sciences
//          Georgia Institute of Technology
//          
//          See AUTHORS file for more information
// 
// Contact: Lavanya Rishishwar
//          lavanya.rishishwar@gatech.edu
//          School of Biological Sciences
//          Georgia Institute of Technology
// ============================================================================ 
// Patent information
// Espitia, H., Chande, A. T., Jordan, I. K., & Rishishwar, L. (2017). 
//      A method of sequence typing with in silico aptamers from a next
//      generation sequencing platform. Google Patents. Patent application 
//      US15/726,005. Retrieved from 
//      https://patents.google.com/patent/US20190108308A1
// ############################################################################


#ifndef _DIR_UTILS_H_
#define _DIR_UTILS_H_

// Source code taken from 
//  - https://insanecoding.blogspot.com/2007/11/pathmax-simply-isnt.html
//  - https://insanecoding.blogspot.com/2007/11/implementing-realpath-in-c.html

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <errno.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

// #include "jlss.h"
// #include "emalloc.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */
// #include "sysstat.h"    /* Fix up for Windows - inc mode_t */

using namespace std;

// =========================================================================
inline bool getCwd(std::string& path)
{
    typedef std::pair<dev_t, ino_t> file_id;

    bool success = false;
    int start_fd = open(".", O_RDONLY); //Keep track of start directory, so can jump back to it later
    if (start_fd != -1)
    {
        struct stat sb;
        if (!fstat(start_fd, &sb))
        {
            file_id current_id(sb.st_dev, sb.st_ino);
            if (!stat("/", &sb)) //Get info for root directory, so we can determine when we hit it
            {
                std::vector<std::string> path_components;
                file_id root_id(sb.st_dev, sb.st_ino);

                while (current_id != root_id) //If they're equal, we've obtained enough info to build the path
                {
                    bool pushed = false;

                    if (!chdir("..")) //Keep recursing towards root each iteration
                    {
                        DIR *dir = opendir(".");
                        if (dir)
                        {
                            dirent *entry;
                            while ((entry = readdir(dir))) //We loop through each entry trying to find where we came from
                            {
                                if ((strcmp(entry->d_name, ".") && strcmp(entry->d_name, "..") && !lstat(entry->d_name, &sb)))
                                {
                                    file_id child_id(sb.st_dev, sb.st_ino);
                                    if (child_id == current_id) //We found where we came from, add its name to the list
                                    {
                                        path_components.push_back(entry->d_name);
                                        pushed = true;
                                        break;
                                    }
                                }
                            }
                            closedir(dir);

                            if (pushed && !stat(".", &sb)) //If we have a reason to contiue, we update the current dir id
                            {
                                current_id = file_id(sb.st_dev, sb.st_ino);
                            }
                        }//Else, Uh oh, can't read information at this level
                    }
                    if (!pushed) { break; } //If we didn't obtain any info this pass, no reason to continue
                }

                if (current_id == root_id) //Unless they're equal, we failed above
                {
                    //Built the path, will always end with a slash
                    path = "/";
                    for (std::vector<std::string>::reverse_iterator i = path_components.rbegin(); i != path_components.rend(); ++i)
                    {
                        path += *i + "/";
                    }
                    success = true;
                }
                fchdir(start_fd);
            }
        }
        close(start_fd);
    }

    return (success);
}
// =========================================================================
inline void relativeDirBaseSplit(std::string& dir, 
                                 std::string& base,
                                 const std::string& path) 
{
  std::string::size_type slash_pos = path.rfind("/"); //Find the last slash
  if (slash_pos != std::string::npos) //If there is a slash
  {
    slash_pos++;
    dir = path.substr(0, slash_pos); //Directory is before slash
    base = path.substr(slash_pos); //And obviously, the file is after
  }
  else //Otherwise, there is no directory present
  {
    dir.clear(); 
    base = path;
  }
}

// =========================================================================
inline bool chdir_getCwd(std::string& path, const std::string& dir)
{
  bool success = false;
  int start_fd = open(".", O_RDONLY); //Open current directory so we can save a handle to it
  if (start_fd != -1)
  {
    if (!chdir(dir.c_str())) //Change to directory
    {
      success = getCwd(path); //And get its path
      fchdir(start_fd); //And change back of course
    }
    close(start_fd);
  }
  return(success);
}

// =========================================================================
static inline bool realpath_file(std::string& resolved_path, const std::string& path)
{
  bool success = false;
  std::string dir;
  std::string base;
  relativeDirBaseSplit(dir, base, path);

  //If there is a directory, get the realpath() for it, otherwise the current directory
  if (dir.size() ? chdir_getCwd(resolved_path, dir) : getCwd(resolved_path))
  {
    resolved_path += base;
    success = true;
  }
  return(success);
}
// =========================================================================
inline bool readlink_internal(std::string& buffer, const std::string& path, ssize_t length)
{
  bool success = false;
  if (length > 0)
  {
    char *buf = new(nothrow) char[length+1]; //Room for Null
    if (buf)
    {
      ssize_t amount = ::readlink(path.c_str(), buf, length+1); //Give room for failure
      if ((amount > 0) && (amount <= length)) //If > length, it was modified mid check
      {
        buf[amount] = 0;
        buffer = buf;
        success = true;
      }
      delete[] buf;
    }
  }
  return(success);
}
// =========================================================================
inline void build_path_base_swap(std::string &path, const std::string& newbase)
{
  string dir;
  string base;
  relativeDirBaseSplit(dir, base, path);

  if (dir.size())
  {
    path = dir + newbase;
  }
  else
  {
    path = newbase;
  }
}
// =========================================================================
inline bool symlink_resolve(std::string& end, const std::string& start)
{
  typedef std::pair<dev_t, ino_t> file_id;

  bool success = false;
  if (start.size())
  {
    std::string path = start; //Need a modifyable copy
    struct stat sb;
    std::set<file_id> seen_links;

    bool resolved_link;
    do //The symlink resolve loop
    {
      resolved_link = false;
      if (!lstat(path.c_str(), &sb))
      {
        file_id current_id(sb.st_dev, sb.st_ino);
        if (seen_links.find(current_id) == seen_links.end()) //Not a link we've seen
        {
          seen_links.insert(current_id); //Add to our set

          if (S_ISLNK(sb.st_mode)) //Another link
          {
            std::string newpath;
            if (readlink_internal(newpath, path, sb.st_size))
            {
              if (newpath[0] == '/') //Absolute
              {
                path = newpath;
              }
              else //We need to calculate the relative path in relation to the current
              {
                build_path_base_swap(path, newpath);
              }
              resolved_link = true;
            } //Else, Link can't be read, time to quit
          }
          else //Yay, it's not a link! got to the last part finally!
          {
            success = realpath_file(end, path);
          }
        } //Else, Nice try, someone linked a link back into a previous link during the scan to try to trick us into an infinite loop
      } //Else, Dangling link, can't resolve
    } while (resolved_link);
  }
  return(success);
}
// =========================================================================
inline bool realpath(std::string& resolved_path, const std::string& path, bool resolve_link = true)
{
  bool success = false;
  if (path.size())
  {
    struct stat sb;
    if (!stat(path.c_str(), &sb))
    {
      bool (*rp)(std::string&, const std::string&) = resolve_link ? symlink_resolve : realpath_file;
      success = S_ISDIR(sb.st_mode) ? chdir_getCwd(resolved_path, path) : rp(resolved_path, path);
    }
  }
  return(success);
}
// =========================================================================
// Source code taken from 
// - https://stackoverflow.com/questions/675039/how-can-i-create-directory-tree-in-c-linux
// - https://github.com/jleffler/soq
/*
@(#)Purpose:        Create all directories in path
@(#)Author:         J Leffler
@(#)Copyright:      (C) JLSS 1990-91,1997-98,2001,2005,2008,2012
*/
inline static int do_mkdir(const char *path, mode_t mode)
{
    struct stat st;
    int         status = 0;

    if (stat(path, &st) != 0)
    {
        /* Directory does not exist. EEXIST for race condition */
        if (mkdir(path, mode) != 0 && errno != EEXIST)
            status = -1;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = -1;
    }

    return(status);
}
// =========================================================================
/**
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards.
*/
inline int mkpath(const char *path, mode_t mode)
{
    char *pp;
    char *sp;
    int   status;
    char *copypath = strdup(path);

    status = 0;
    pp     = copypath;
    while (status == 0 && (sp = strchr(pp, '/')) != 0)
    {
        if (sp != pp)
        {
            /* Neither root nor double slash in path */
            *sp    = '\0';
            status = do_mkdir(copypath, mode);
            *sp    = '/';
        }
        pp = sp + 1;
    }
    if (status == 0)
        status = do_mkdir(path, mode);
    free(copypath);
    return (status);
}

#endif