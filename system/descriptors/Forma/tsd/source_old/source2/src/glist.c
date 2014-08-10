/* Insert */
#include "glist.h"

GList* CreateGList(void (*freeF)(void**))
{
    GList* s = (GList*) malloc(sizeof(GList));

    s->freeFunction = freeF;
    s->first = s->last = NULL;

    return s;
}

void GListInsert (GList * qp, void* t)
{
    if(qp == NULL) return;
   GListNode * n = (GListNode *) calloc(1,sizeof(GListNode));

   /* Check if malloc succeeded */

   if (n == NULL) {
      fprintf(stderr, "Out of memory\n");
      exit(1);
   }

   /* Copy the data */

   n->data = t;
   n->next = NULL;

   /* If the GList was empty, just add this one element */

   if (qp->last == NULL)
   {
      qp->first = qp->last = n;
   }else {
      qp->last->next = n;
      qp->last = n;
   }
}

void DestroyGList(GList** s)
{
    if(s == NULL || *s == NULL) return;

    GListNode* t = (*s)->first;
    int i = 0;
    while( t != NULL)
    {
        if(t->data != NULL)
            (*s)->freeFunction(&(t->data));
        GListNode* temp = t;
        t = t->next;
        free(temp);
        fprintf(stderr,"cleaning up %d\n",i++);
    }

    free(*s);

    *s = NULL;
}
