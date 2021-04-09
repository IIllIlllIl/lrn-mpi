#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{
    int i, j, n;
    scanf("%d", &n);
    char *visit = (char *)malloc(n);
    int *prm = (int *)malloc(n * sizeof(int));

    for (i = 0; i < n; i++)
        visit[i] = 0;
    for (i = 0; i < n; i++)
        prm[i] = 0;

    for (i = 2; i <= n; i++)
    {
        if (!visit[i])
        {
            prm[++prm[0]] = i;
        }
        for (j = 1; j <= prm[0] && i * prm[j] <= n; j++)
        {
            visit[i * prm[j]] = 1;
            if (i % prm[j] == 0)
            {
                break;
            }
        }
    }

    printf("%d\n", *prm);

    int *index;
    index = ++prm;
    for (i = 0; i < 10; i++)
        printf("%d ", *(++index));
}