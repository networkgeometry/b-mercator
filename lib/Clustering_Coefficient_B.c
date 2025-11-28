//This function computes the clustering coefficient of nodes (Type1) and features (Type2) and save them in CCO_n and CCo_f respectively 
//Edeglist is network edgelist
//N_n is the number of nodes
//N_f is the number of features
//N_Links is the number of links 


void Clustering_Coefficient_B(int **Edgelist, int N_n, int N_f,  int N_Links, float * CCO_n , float * CCO_f)
{

    double Clustering_Coeff_n[N_n][2];
    double Clustering_Coeff_f[N_f][2];

    struct Nodes {
        int Label;
        struct Nodes* next;
    };

    struct Nodes* NL_N[N_n]; //if the network is too big then its better to use malloc to allocate the memory for NL_N and NL_F
    struct Nodes* NL_F[N_f]; //NL_N contains neighbors lists for nodes and NL_F contains neighbors list for features 
    struct Nodes* temp;
    for(int i=0;i<N_n;i++)   // make a linked-list of the neighbors of Nodes
    {
        NL_N[i]=0;
    }
    

    for (int i=0; i<N_n ; i++) //For node of Type1 i
    {
        for(int j=0; j<N_Links; j++) //For each link j
        {
            if(Edgelist[j][0]==i) //if it is a link like i---x
            {
                //add x to the neighbors of i in the linked list
                temp = (struct Nodes*)malloc(sizeof(struct Nodes));
                temp->Label =Edgelist[j][1];

                if(NL_N[i]==0) //if it is the first neighbors added to the linked list
                 {
                     NL_N[i]=temp;
                     temp->next = 0;
                 }
                 else
                 {
                     temp->next= NL_N[i];
                     NL_N[i]= temp;
                 }
            }
        }
    }

   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~// Make a linked-list of the neighbors of features
    for(int i=0;i<N_f;i++)
    {
        NL_F[i]=0;
    }
    for (int i= N_n; i<N_n + N_f ; i++) //the label of features are in [N_n , N_n + N_f]
    {
        //printf("%d\n", i);
        for(int j=0; j< N_Links; j++)
        {
            if(Edgelist[j][1]==i)
            {
                temp = (struct Nodes*)malloc(sizeof(struct Nodes));
                temp->Label =Edgelist[j][0];

                if(NL_F[i-N_n]==0) //The first element
                 {

                     NL_F[i-N_n]=temp;
                     temp->next = 0;
                 }
                 else
                 {

                     temp->next= NL_F[i-N_n];
                     NL_F[i-N_n]= temp;
                 }
            }
        }
    }
   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    struct Nodes * temp3;
    for (int i=0; i<N_n ; i++)//For every node of Type 1 called i do
    {
        int num_tri=0;  //it keeps the number of triangels
        int count_Ni=0;
        temp3= NL_N[i];

        while(temp3 != 0)  //The number of neighbores of node i
        {
            count_Ni= count_Ni+1;
            temp3= temp3->next;
        }

        int N_i[count_Ni];
        temp3= NL_N[i];


        for (int cc=0; cc<count_Ni ; cc++) //The set of neighbors of i
        {
            N_i[cc]= temp3->Label;
            temp3= temp3->next;

        }


        for (int l=0; l<(count_Ni-1); l++) //For each neighbor of i called l
        {

            int ind_l= N_i[l]-N_n;

            int count_Nl=0;
            struct Nodes * temp4;
            temp4= NL_F[ind_l];
            while(temp4 != 0)  //The number of neighbores of node l
            {
                count_Nl= count_Nl+1;
                temp4= temp4->next;
            }

            int N_l[count_Nl];
            temp4= NL_F[ind_l];

            for (int cc=0; cc<count_Nl ; cc++) //The set of neighbors of node l
            {
                N_l[cc]= temp4->Label;
                temp4= temp4->next;

            }

            for (int m=(l+1); m < count_Ni; m++) //for every  neighbor of i called m do
            {
                int ind_m= N_i[m]-N_n;
                int count_Nm=0;
                temp4= NL_F[ind_m];
                while(temp4 != 0)  //The number of neighbores of node m
                {
                    count_Nm= count_Nm+1;
                    temp4= temp4->next;
                }
                int N_m[count_Nm];
                temp4= NL_F[ind_m];

                for (int cc=0; cc<count_Nm ; cc++) //the set of neighbors of m
                {
                    N_m[cc]= temp4->Label;
                    temp4= temp4->next;

                }

                int Num_CN= 0;
                for (int ii=0; ii<count_Nl; ii++)  //the number of common neighbors between l and m
                {
                    for (int jj=0; jj < count_Nm ; jj++)
                    {
                        if (N_l[ii]==N_m[jj])
                        {
                            Num_CN=Num_CN + 1;
                        }
                    }

                }

                if(Num_CN >1) // it l and m share common neighbors except i
                {
                   num_tri= num_tri + 1; //icrease the number of triangels going through i by 1
                }

            }
        }

       Clustering_Coeff_n[i][0]=i;

       double cco=0;  // cco = 2T/ (ki * (ki -1))
       double A= 2*num_tri;
       double B = count_Ni * (count_Ni - 1);

       if(B>0)
       {
           cco= A/B;
       }
       else
       {
           cco=0;
       }
       Clustering_Coeff_n[i][1]= cco;

    }
    //~~~~~~~~~~~~~~~~~~~Clustering Coefficient of nodes of Type2
    for (int j=0; j<N_f ; j++)//For every node of Type 2 do
    {

        int num_tri=0;  //it keeps the number of triangels
        int count_Nj=0;
        temp3= NL_F[j];

        while(temp3 != 0)  //The number of neighbores of node N_i
        {
            count_Nj= count_Nj+1;
            temp3= temp3->next;
        }

        int N_j[count_Nj];
        temp3= NL_F[j];

        for (int cc=0; cc<count_Nj ; cc++)
        {
            N_j[cc]= temp3->Label;
            temp3= temp3->next;
        }

        for (int l=0; l<(count_Nj-1); l++)
        {

            int ind_l= N_j[l];
            int count_Nl=0;
            struct Nodes * temp4;
            temp4= NL_N[ind_l];
            while(temp4 != 0)  //The number of neighbores of node N_l
            {
                count_Nl= count_Nl+1;
                temp4= temp4->next;
            }

            int N_l[count_Nl];
            temp4= NL_N[ind_l];

            for (int cc=0; cc<count_Nl ; cc++)
            {
                N_l[cc]= temp4->Label;
                temp4= temp4->next;

            }

            for (int m=(l+1); m < count_Nj; m++)
            {
                int ind_m= N_j[m];
                int count_Nm=0;
                temp4= NL_N[ind_m];
                while(temp4 != 0)  //The number of neighbores of node N_m
                {
                    count_Nm= count_Nm+1;
                    temp4= temp4->next;
                }
                int N_m[count_Nm];
                temp4= NL_N[ind_m];
                for (int cc=0; cc<count_Nm ; cc++)
                {
                    N_m[cc]= temp4->Label;
                    temp4= temp4->next;
                }

                int Num_CN= 0;
                for (int ii=0; ii<count_Nl; ii++)
                {
                    for (int jj=0; jj < count_Nm ; jj++)
                    {
                        if (N_l[ii]==N_m[jj])
                        {
                            Num_CN=Num_CN + 1;
                        }
                    }

                }

                if(Num_CN >1)
                {
                    num_tri= num_tri + 1;
                }

            }
        }


        Clustering_Coeff_f[j][0]=j+N_n;
        double cco=0;
        double A= 2*num_tri;
        double B = count_Nj * (count_Nj - 1);
        if (B>0)
        {
            cco= A/B;
        }
        else
        {
            cco=0;
        }
        Clustering_Coeff_f[j][1]= cco;
    }
    //~~~~~~~~~~~~~~~~~~~~~~~~~
    float Avg_CCO_N=0;
    for(int a=0; a<N_n; a++)
    {
        Avg_CCO_N= Avg_CCO_N + Clustering_Coeff_n[a][1];
    }
    * CCO_n = Avg_CCO_N / N_n ;
    //****************************Compute the average CCO of Nodes of Type 2 for the network
    float Avg_CCO_F=0;
    for(int a=0; a<N_f; a++)
    {
        Avg_CCO_F= Avg_CCO_F + Clustering_Coeff_f[a][1];
    }

    * CCO_f = Avg_CCO_F / N_f ;


}

