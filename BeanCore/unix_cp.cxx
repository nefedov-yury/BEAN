#if defined (__unix__)

// copy whole directory
// Use to copy sqlite data files to shared memory or /tmp

#include <cstdio>
#include <cstdlib>
#include <cerrno>

// opendir,mkdir,closedir,remove
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

// default permissions for new directories and files:
#define DIR_MODE        (FILE_MODE | S_IXUSR | S_IXGRP | S_IXOTH)
#define FILE_MODE       (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)

#define BUFSIZE 4096

using namespace std;

//--------------------------------------------------------------------
static void err_prt(const char *s1, const char *s2)
//--------------------------------------------------------------------
{
  fprintf(stderr,"unix_cp error: %s ", s1);
  perror(s2);
  exit(1);
}

//--------------------------------------------------------------------
void copy_file(const char *file_in, const char* file_out)
//--------------------------------------------------------------------
{
   char buf[BUFSIZE];

   // open file
   FILE* fp_in = fopen(file_in,"r");
   if( !fp_in ) err_prt("can not open ", file_in);

   // create new file
   FILE* fp_out = fopen(file_out,"w");
   if( !fp_out ) err_prt("can not create ", file_out);

   // copy
   size_t nchars;
   while( (nchars = fread(buf,1, BUFSIZE,fp_in)) > 0 ) {
        if( nchars != fwrite(buf,1,nchars,fp_out) ) {
          err_prt("write to ", file_out);
        }
        if( nchars != BUFSIZE ) {
          if( feof(fp_in) ) break;
          err_prt("read from ", file_in);
        }
   }

   // close files
   if( fclose(fp_in) ) err_prt("can not close ", file_in);
   if( fclose(fp_out) ) err_prt("can not close ", file_out);
}

//--------------------------------------------------------------------
void copy_dir(const char *dir_in, const char* dir_out)
//--------------------------------------------------------------------
{
   char old_file[257];
   char new_file[257];

   // open dir
   DIR* dp_in = opendir(dir_in);
   if( !dp_in ) err_prt("can not open dir ", dir_in);

   // create new dir -> mkdtemp do it
//    int dp_out = mkdir(dir_out,DIR_MODE);
//    if( dp_out < 0 ) err_prt("can not create dir ", dir_out);

   // copy regular files
   struct dirent* de;
   while( (de = readdir(dp_in)) != NULL ) {
      // get file information
      snprintf(old_file,sizeof(old_file),"%s/%s",dir_in,de->d_name);
      struct stat statbuf;
      if( stat(old_file,&statbuf) == -1 ) continue;

      if( S_ISREG(statbuf.st_mode) ) { // man stat.h
        snprintf(new_file,sizeof(new_file),"%s/%s",dir_out,de->d_name);
        // printf(" %s -> %s\n",old_file, new_file);
        copy_file(old_file, new_file);
      }
   }

   // close dir
   if( closedir(dp_in) < 0 ) err_prt("can not close ",dir_in);
}

//--------------------------------------------------------------------
void rm_whole_dir(const char *dir)
//--------------------------------------------------------------------
{
   char file[257];

   // open dir
   DIR* dp = opendir(dir);
   if( !dp ) err_prt("can not open dir ", dir);

   // remove files
   struct dirent* de;
   while( (de = readdir(dp)) != NULL ) {
      snprintf(file,sizeof(file),"%s/%s",dir,de->d_name);
      struct stat statbuf;
      if( stat(file,&statbuf) == -1 ) continue;

      // assume that all files are regular
      if( S_ISREG(statbuf.st_mode) ) { // man stat.h
        if( remove(file) ) err_prt("can not remove ",file);
      }
   }

   // close dir and remove it
   if( closedir(dp) < 0 ) err_prt("can not close ",dir);
   if( remove(dir) ) err_prt("can not remove ",dir);
}

//--------------------------------------------------------------------
const char* copy_dir_temp(const char *dir_in, const char* basedir)
//--------------------------------------------------------------------
{
   static char tmp_dir[256];
   snprintf(tmp_dir,sizeof(tmp_dir),"%s/bean_XXXXXX",basedir);
//    mktemp(tmp_dir); unsafe -> use mkdtemp
   if ( mkdtemp(tmp_dir) == NULL ) {
      err_prt("can not create dir ",tmp_dir);
   }
   copy_dir(dir_in, tmp_dir);

   return tmp_dir;
}

#endif
