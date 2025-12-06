    //FOR DEBUG
      /*if(ipath==158){
            printf("path ---> (");
            for(int j=0; j<length; j++){
                printf("%d ", this_path->dir[j]);
            }
            printf(")\n");
        int countertemp = 0;
        for (size_t i = 0; i < 10; ++i) {
          countertemp++;
          int test = i * 12345 % SITES_ON_NODE; // pseudo-random
          //int nbrs = walk_dir(test, dir);
          int x,y,z,t;
          coords_from_site_index(test, &x,&y,&z,&t);
          //coords_from_site_index(nbrs,&xn,&yn,&zn,&tn);
          printf("i=%d site=(%d,%d,%d,%d) "
              "\n",
              test, x,y,z,t);

          printf("oprod_along_path[%d][%d] = ", length, test);print_su3(&oprod_along_path[length][test]);
          if(countertemp==10){
              printf("countertemp = %d\n", countertemp); 
          }
        }
        exit(0);
      }*/
     /*if(ipath == 98){
            int countertemp = 0;
            printf("path ---> (");
            for(int j = 0; j < length; j++){
                printf("%d ", this_path->dir[j]);
            }
            printf(")\n");

            for(size_t i = 0; i < 10; ++i){
                countertemp++;
                int test = i * 12345 % SITES_ON_NODE;
                int x, y, z, t;
                coords_from_site_index(test, &x, &y, &z, &t);

                printf("i=%d site=(%d,%d,%d,%d)\n", test, x, y, z, t);

                printf("mats_along_path[%d][%d] = ", length, test);
                print_su3(&mats_along_path[length][test]);

                if(countertemp == 10){
                    printf("countertemp = %d\n", countertemp);
                }
            }
            exit(0);
        }*/

        /*if(ipath == 99){
            int countertemp = 0;
            printf("path ---> (");
            for(int j = 0; j < length; j++){
                printf("%d ", this_path->dir[j]);
            }
            printf(")\n");

            for(int d = 0; d < 4; d++){
                printf("force_accum direction %d:\n", d);
                countertemp = 0;
                for(size_t i = 0; i < 10; ++i){
                    countertemp++;
                    int test = i * 12345 % SITES_ON_NODE;
                    int x, y, z, t;
                    coords_from_site_index(test, &x, &y, &z, &t);

                    printf("i=%d site=(%d,%d,%d,%d)\n", test, x, y, z, t);
                    printf("force_accum[%d][%d] = ", d, test);
                    print_su3(&force_accum[d][test]);

                    if(countertemp == 10){
                        printf("countertemp = %d\n", countertemp);
                    }
                }
                printf("\n");
            }
            exit(0);
        }*/
        //DEBUG
    
    /*int countertemp = 0;

    for (int dir = 0; dir < 4; dir++) {
        printf("mom direction %d:\n", dir);
        countertemp = 0;

        for (size_t n = 0; n < 10; ++n) {
            countertemp++;
            int test = (int)(n * 12345 % SITES_ON_NODE);

            int x, y, z, t;
            coords_from_site_index(test, &x, &y, &z, &t);

            printf("i=%d site=(%d,%d,%d,%d)\n", test, x, y, z, t);

            anti_hermitmat *M = &mom[test][dir];

            printf("mom[%d][%d] = ", dir, test);
            printf("(%g,%g) (%g,%g) (%g,%g)  %g %g %g\n",
                   M->m01.real, M->m01.imag,
                   M->m02.real, M->m02.imag,
                   M->m12.real, M->m12.imag,
                   M->m00im, M->m11im, M->m22im);

            if (countertemp == 10) {
                printf("countertemp = %d\n", countertemp);
            }
        }
        printf("\n");
    }

    exit(0);*/