#ifndef AMTOOLBOX_H
#define AMTOOLBOX_H

typedef struct adaptloopstatevar
{
  int    loops;
  int    nsigs;
  double *state;
  double corr;
  double mult;
  double limit;
  double minlvl;
  double *a1;
  double *factor, *expfac, *offset;
} adaptloopstate;


#endif /* AMTOOLBOX_H */
