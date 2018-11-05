/*
 * GET SIGNAL VALUE, (original 2016/06/30, R.Yoshida)
 * setting IFAIL to 0 to indicate not appropriate value.
 */

#include <signal.h>

int get_sigint(int *ifail)
{
  if (SIGINT < 0) {
      *ifail = 1;
      return 0;
    }
  else {
      *ifail = 0;
      return SIGINT;
    }
}

int get_sigquit(int *ifail)
{
  if (SIGQUIT < 0) {
      *ifail = 1;
      return 0;
    }
  else {
      *ifail = 0;
      return SIGQUIT;
    }
}

int get_sigabrt(int *ifail)
{
  if (SIGABRT < 0) {
      *ifail = 1;
      return 0;
    }
  else {
      *ifail = 0;
      return SIGABRT;
    }
}

int get_sigfpe(int *ifail)
{
  if (SIGFPE < 0) {
      *ifail = 1;
      return 0;
    }
  else {
      *ifail = 0;
      return SIGFPE;
    }
}

int get_sigsegv(int *ifail)
{
  if (SIGSEGV < 0) {
      *ifail = 1;
      return 0;
    }
  else {
      *ifail = 0;
      return SIGSEGV;
    }
}

int get_sigterm(int *ifail)
{
  if (SIGTERM < 0) {
      *ifail = 1;
      return 0;
    }
  else {
      *ifail = 0;
      return SIGTERM;
    }
}

