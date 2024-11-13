
#include<iostream>
#include "Algorithm/zoxc.h"
#include "Algorithm/zgrouping.h"
extern "C"{
#include "Jerasure-1.2A/jerasure.h"
#include "Jerasure-1.2A/cauchy.h"
}
#include<ctime>
#include<unistd.h>
#include <sys/time.h>
#include<assert.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<bitset>
#include<vector>
#include<set>
#include<algorithm>
using namespace std;

long long T0 = 0;
long long T1 = 0;

long long COM_XOR_TIME = 0;
long long COM_CPY_TIME = 0;

// 统计一个矩阵中一的个??
int get_sum_ones(int* matrix, int k, int m, int w) {
    int sum = 0;
    for(int i = 0; i < k * m; i++) {
        sum += cauchy_n_ones(matrix[i], w);
    }
    return sum;
}
// 统计智能调度的xor+copy的总次??
int get_sum_schedules(int** schedule) {
    int i = 0;
    while(schedule[i][0] >= 0) {
        i++;
    }
    return i;
}
// 统计公共异或操作中xor的总次??
int get_sum_common(vector<int*>& schedule) {
    int cnt = 0;
    for(int i = 0; i < schedule.size(); i++) {
        if(schedule[i][2] % 2 == 1) {
            cnt++;
        }
    }
    return cnt;
}
// 输出智能调度序列
void print_smart_schedule(int** schedule) {
    int i = 0;
    while(schedule[i][0] >= 0) {
        for(int j = 0; j < 5; j++) {
            cout << schedule[i][j] << " ";
        }
        cout << endl;
        i++;
    }
}
// 输出公共异或调度序列
void print_comxor_schedule(vector<int*> schedule) {
    for(int i = 0; i < schedule.size(); i++) {
        for(int j = 0; j < 3; j++) {
            cout << schedule[i][j] << " ";
        }
        cout << endl;
    }
}
// 对位矩阵进行划分
vector<int*> split(int* bitmatrix, int k, int m, int w, int nr) {
    vector<int*> ans(nr);
    int new_k = k / nr;
    for(int i = 0; i < nr; i++) {
        ans[i] = new int[new_k * m * w * w];
    }
    for(int i = 0; i < k * m * w * w; i++) {
        int bit_x = i / (k * w);
        int bit_y = i % (k * w);
        int ans_i = bit_y / (new_k * w);
        int ans_j = bit_y % (new_k * w) + bit_x * new_k * w;
        ans[ans_i][ans_j] = bitmatrix[i];
    }
    return ans;
}

// 统计智能调度中异或的次数
int get_xor_schedules(int** schedule) {
    int i = 0;
    int Xor = 0;
    while(schedule[i][0] >= 0) {
        if(schedule[i][4] == 1) {
            Xor++;
        }
        i++;
    }
    return Xor;
}
// 统计智能调度中copy的次??
int get_copy_schedules(int** schedule) {
    int i = 0;
    int copy = 0;
    while(schedule[i][0] >= 0) {
        if(schedule[i][4] == 0) {
            copy++;
        }
        i++;
    }
    return copy;
}

// 统计公共异或操作中xor的总次??
int get_xor_common(vector<int*>& schedule) {
    int cnt = 0;
    for(int i = 0; i < schedule.size(); i++) {
        if(schedule[i][2] % 2 == 1) {
            cnt++;
        }
    }
    return cnt;
}
// 统计公共异或操作中的copy次数
int get_copy_common(vector<int*>& schedule) {
    int cnt = 0;
    for(int i = 0; i < schedule.size(); i++) {
        if(schedule[i][2] % 2 == 0) {
            cnt++;
        }
    }
    return cnt;
}


// 测试智能调度和公共异或，对于划分后的效果————时间和异或次数的减??
void test1() {
    const int k = 24;
    const int m = 3;
    const int w = 8;
    const int nr = 1;
    const int new_k = k / nr;
    printf("**************k = %d, m = %d, w = %d, nr = %d**************\n", k, m, w, nr);
    
    int* cauchy = cauchy_good_general_coding_matrix(k, m, w);
    jerasure_print_matrix(cauchy, m, k, w);
    int ones = get_sum_ones(cauchy, k, m , w);
    cout << "bitmatrix ones: \t\t" << ones << endl;
    cout << "xor nums: \t\t\t" << ones - m * w << endl;

    int* bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, cauchy);
    // jerasure_print_bitmatrix(bitmatrix, m * w, k * w, w);


    vector<int*> split_bitmatrix = split(bitmatrix, k, m, w, nr);
    // for(int i = 0; i < split_bitmatrix.size(); i++) {
    //     cout << "**********************" << endl;
    //     jerasure_print_bitmatrix(split_bitmatrix[i], m * w , k / nr * w, w);
    // }


    

    // int** schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
    // cout << "schedlue xor nums:\t\t" << get_sum_schedules(schedule) - m * w << endl;
    // print_smart_schedule(schedule);

    int xor_num = 0;
    for(int i = 0; i < nr; i++) {
        int** schedule = jerasure_smart_bitmatrix_to_schedule(new_k, m, w, split_bitmatrix[i]);
        xor_num += get_sum_schedules(schedule) - m * w;
        // print_smart_schedule(schedule);
        // cout << "**************************" << endl;
        free(schedule);
    }
    cout << "split---schedlue xor nums:\t" << xor_num + (nr - 1) * m * w << endl;


    // ZOXC xoc(k * w, m * w);
    // clock_t start1 = clock();
    // int n_XOR = xoc.grouping_1s(bitmatrix, false);
    // clock_t end1 = clock();
    // cout << "common xor nums: \t\t" << get_sum_common(xoc.schedule) << "---------->" << (double)(end1 - start1) / CLOCKS_PER_SEC << "s" << endl;
    // cout << "n_XOR: \t\t\t\t" << n_XOR << endl;
    
    // print_comxor_schedule(xoc.schedule);
    // cout << "************************************" << endl;
    // print_comxor_schedule(xoc.intermedia_schedule);
    // cout << xoc.schedule.size() << endl;

    int xor_num_split = 0;
    int n_xor_split = 0;
    clock_t start2 = clock();
    for(int i = 0; i < nr; i++) {
        ZOXC xoc_split(new_k * w, m * w);
        n_xor_split += xoc_split.grouping_1s(split_bitmatrix[i], false);
        xor_num_split += get_sum_common(xoc_split.schedule);
        // print_comxor_schedule(xoc_split.intermedia_schedule);
        // cout << "**************************" << endl;
    }
    clock_t end2 = clock();
    cout << "split---common xor nums: \t" << xor_num_split + (nr - 1) * m * w << "---------->" << (double)(end2 - start2) / CLOCKS_PER_SEC << "s" << endl;
    // cout << "split----n_XOR: \t\t\t\t" << n_XOR << endl;

}

void schedule_encode_single_chunk(char* data, int len, char**& parities, int K, int M, int W, int packetsize, int** en_schedule) {
    assert(data != NULL);
    assert(parities != NULL);
    assert(len ==  packetsize * K * W);
    
    // struct timeval t0,t1;
    struct timespec t0, t1;
    t0.tv_nsec = 0;
    t0.tv_sec = 0;
    t1.tv_nsec = 0;
    t1.tv_sec = 0;
    
    char **tdata = (char**)malloc(K*sizeof(char*));
    for(int i = 0;i<K;i++)
        tdata[i] = data + packetsize * i * W;
    // gettimeofday(&t0,NULL);
    clock_gettime(CLOCK_REALTIME, &t0);
    jerasure_schedule_encode(K,M,W,en_schedule,tdata,parities, W*packetsize, packetsize);
    clock_gettime(CLOCK_REALTIME, &t1);
    // gettimeofday(&t1,NULL);
    // T0 += diff_us(t0, t1);
    // cout << "T0-------------------------" << T0 << endl;
    T0 += (t1.tv_sec-t0.tv_sec)*1000000000+(t1.tv_nsec-t0.tv_nsec);
    free(tdata);
}

void do_scheduled_operations(vector<int*> &schedule, vector<char*>& intermedia, char **&data, char **&parities, int K, int M, int W, int packetsize)
{
    char *sptr=NULL;
    char *dptr=NULL;
    int sch_len = schedule.size();
    int s,d,op,d_idx,d_off, s_idx, s_off;
    struct timespec t0, t1;
    t0.tv_nsec = 0;
    t0.tv_sec = 0;
    t1.tv_nsec = 0;
    t1.tv_sec = 0;
    for(int i = 0;i<sch_len;i++)
    {
        s = schedule[i][0];
        d = schedule[i][1];
        op = schedule[i][2];
        s_idx = s / W;
        s_off = (s % W)*packetsize;
        d_idx = d / W;
        d_off = d % W * packetsize;
        switch(op)
        {
        case 0:
            sptr = data[s_idx] + s_off;
            dptr = parities[d_idx] + d_off;
            clock_gettime(CLOCK_REALTIME, &t0);
            memcpy(dptr, sptr, packetsize);
            clock_gettime(CLOCK_REALTIME, &t1);
            COM_CPY_TIME += (t1.tv_sec-t0.tv_sec)*1000000000+(t1.tv_nsec-t0.tv_nsec);
            break;
        case 1:
            sptr = data[s_idx] + s_off;
            dptr = parities[d_idx] + d_off;
            clock_gettime(CLOCK_REALTIME, &t0);
            fast_xor(sptr, dptr, dptr, packetsize);
            clock_gettime(CLOCK_REALTIME, &t1);
            COM_XOR_TIME += (t1.tv_sec-t0.tv_sec)*1000000000+(t1.tv_nsec-t0.tv_nsec);
            break;
        case 2:
            sptr = intermedia[s];
            dptr = parities[d_idx] + d_off;
            clock_gettime(CLOCK_REALTIME, &t0);
            memcpy(dptr, sptr, packetsize);
            clock_gettime(CLOCK_REALTIME, &t1);
            COM_CPY_TIME += (t1.tv_sec-t0.tv_sec)*1000000000+(t1.tv_nsec-t0.tv_nsec);
            break;
        case 3:
            sptr = intermedia[s];
            dptr = parities[d_idx] + d_off;
            clock_gettime(CLOCK_REALTIME, &t0);
            fast_xor(sptr, dptr, dptr, packetsize);
            clock_gettime(CLOCK_REALTIME, &t1);
            COM_XOR_TIME += (t1.tv_sec-t0.tv_sec)*1000000000+(t1.tv_nsec-t0.tv_nsec);
            break;
        case 4:
            sptr = data[s_idx] + s_off;
            dptr = intermedia[d];
            clock_gettime(CLOCK_REALTIME, &t0);
            memcpy(dptr, sptr, packetsize);
            clock_gettime(CLOCK_REALTIME, &t1);
            COM_CPY_TIME += (t1.tv_sec-t0.tv_sec)*1000000000+(t1.tv_nsec-t0.tv_nsec);
            break;
        case 5:
            sptr = data[s_idx] + s_off;
            dptr = intermedia[d];
            clock_gettime(CLOCK_REALTIME, &t0);
            fast_xor(sptr, dptr, dptr, packetsize);
            clock_gettime(CLOCK_REALTIME, &t1);
            COM_XOR_TIME += (t1.tv_sec-t0.tv_sec)*1000000000+(t1.tv_nsec-t0.tv_nsec);
            break;
        }
    }
}

void commonxor_encode_single_chunk(char *data, int len, char **&parities, int K, int M, int W, int packetsize, vector<int*>& en_schedule, vector<char*>& intermedia)
{
    assert(data != NULL);
    assert(parities != NULL);
    assert(len == packetsize * K * W);
    struct timespec t0, t1;
    t0.tv_nsec = 0;
    t0.tv_sec = 0;
    t1.tv_nsec = 0;
    t1.tv_sec = 0;
    char **tdata = (char**)malloc(K*sizeof(char*));
    for(int i = 0;i<K;i++)
        tdata[i] = data + packetsize * i * W;
    clock_gettime(CLOCK_REALTIME, &t0);
    do_scheduled_operations(en_schedule, intermedia, tdata, parities, K, M, W, packetsize);
    clock_gettime(CLOCK_REALTIME, &t1);
    T1 += (t1.tv_sec-t0.tv_sec)*1000000000+(t1.tv_nsec-t0.tv_nsec);
    // cout << "T1:------------------------------ " << T1 << endl; 
    free(tdata);
}

// 测试智能调度和公共异或，编码吞吐量的对比,10MB的数??
void test2(const int k, int m, const int w, const int nr , const int packetsize, long long data_size) {
    const int new_k = k / nr;
    const int blocksize = packetsize * k  *w;
    struct timeval t0,t1;
    long long  diff;
    while(data_size % blocksize != 0) {
        data_size++;
    }
    clock_t start, end;
    printf("**************k = %d, m = %d, w = %d, nr = %d, size = %d, packetsize = %d**************\n", k, m, w, nr, data_size, packetsize);

    int* cauchy = cauchy_good_general_coding_matrix(k, m, w);           // 生成编码矩阵
    // cauchy = cauchy + k;
    // m = m - 1;
//     int cauchy[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
// 66, 235, 38, 13, 138, 73, 1, 147, 103, 58, 93, 161, 163, 24, 42, 247, 67, 117, 181, 60, 234, 175, 238, 15, 151, 108, 85, 18, 80, 190, 88, 109, 118, 71, 229, 167, 182, 6, 206, 29, 57, 172, 27, 102, 156, 107, 2, 199, 
// 12, 183, 15, 14, 8, 37, 134, 97, 218, 253, 9, 234, 177, 158, 76, 172, 147, 161, 11, 53, 122, 216, 133, 52, 255, 56, 98, 228, 89, 42, 41, 239, 64, 80, 217, 70, 155, 235, 132, 236, 209, 1, 200, 68, 21, 104, 195, 108, 
// 179, 79, 254, 40, 184, 48, 100, 91, 95, 152, 166, 230, 68, 125, 203, 236, 156, 19, 89, 87, 41, 28, 143, 60, 14, 238, 150, 3, 94, 164, 15, 57, 140, 109, 163, 255, 112, 237, 104, 234, 24, 1, 18, 200, 76, 134, 51, 81};
    jerasure_print_matrix(cauchy, m, k, w);

    int ones = get_sum_ones(cauchy, k, m , w);                           // 统计矩阵??1的个??
    cout << "bitmatrix ones: \t\t" << ones << endl;
    cout << "xor nums: \t\t\t" << ones - m * w << endl;

    int* bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, cauchy);     // 转换为位矩阵
    jerasure_print_bitmatrix(bitmatrix, m*w, k*w, w);
    
    vector<int*> split_bitmatrix = split(bitmatrix, k, m, w, nr);       // 将位矩阵进行分割
    
    
    // int xor_num = 0;                                                  // 统计使用智能调度后的异或总次??
    vector<int**> split_schedule;               // 保存智能调度，划分后的智能调??                     
    for(int i = 0; i < nr; i++) {
        int** schedule = jerasure_smart_bitmatrix_to_schedule(new_k, m, w, split_bitmatrix[i]);
        // xor_num += get_sum_schedules(schedule) - m * w;
        split_schedule.push_back(schedule);
    }
    // cout << "split---schedlue xor nums:\t" << xor_num + (nr - 1) * m * w << endl;


    // int xor_num_split = 0;
    vector<vector<int*> > xoc_split_schedule;       // 保存划分后的公共异或调度
    vector<int> intermedia_size;
    gettimeofday(&t0,NULL);
    for(int i = 0; i < nr; i++) {
        ZOXC xoc_split(new_k * w, m * w);
        xoc_split.grouping_1s(split_bitmatrix[i], false);
        // xor_num_split += get_sum_common(xoc_split.schedule);
        vector<int*> xoc_schedule(xoc_split.schedule.size());
        for(int i = 0; i < xoc_split.schedule.size(); i++) {
            xoc_schedule[i] = new int[3];
            xoc_schedule[i][0] = xoc_split.schedule[i][0];
            xoc_schedule[i][1] = xoc_split.schedule[i][1];
            xoc_schedule[i][2] = xoc_split.schedule[i][2];
        }
        xoc_split_schedule.push_back(xoc_schedule);
        intermedia_size.push_back(xoc_split.intermedia_schedule.size());
        cout << xoc_split.intermedia_schedule.size() << "size------------------------------------" << endl;
    }
    gettimeofday(&t1,NULL);
    // cout << "split---common xor nums: \t" << xor_num_split + (nr - 1) * m * w << "---------->" << diff_us(t0,t1) << "us" << endl;

    // for(int i = 0; i < split_schedule.size(); i++) {
    //     print_smart_schedule(split_schedule[i]);
    //     printf("***************\n");
    // }

    // for(int i = 0; i < xoc_split_schedule.size(); i++) {
    //     print_comxor_schedule(xoc_split_schedule[i]);
    //     printf("-------------\n");
    // }

    // 恢复调度编号
    int ii = 0;
    while(split_schedule[0][ii][0] != -1) {
        split_schedule[0][ii][2] += (k - new_k);
        if(split_schedule[0][ii][0] >= new_k) {
            split_schedule[0][ii][0] += (k - new_k);
        }
        ii++;
    } 
    for(int i = 1; i < split_schedule.size(); i++) {
        int ii = 0;
        while(split_schedule[i][ii][0] != -1) {
            split_schedule[i][ii][2] += (k - new_k);
            split_schedule[i][ii][4] = 1;
            if(split_schedule[i][ii][0] >= new_k) {
                split_schedule[i][ii][0] += (k - new_k);
            }
            else {
                split_schedule[i][ii][0] = split_schedule[i][ii][0] + i * new_k; 
            }
            ii++;
        }
    }
    
    // for(int i = 0; i < split_schedule.size(); i++) {
    //     print_smart_schedule(split_schedule[i]);
    //     printf("***************\n");
    // }

    // 恢复公共异或的编??
    //      0 for copy data,
    //      1 for xor data,
    //      2 for copy from intermedia,
    //      3 for xor from intermedia,
    //      4 for copy to intermedia,
    //      5 for xor to intermedia
    for(int i = 1; i < xoc_split_schedule.size(); i++) {
        
        for(int j = 0; j < xoc_split_schedule[i].size(); j++) {
            switch (xoc_split_schedule[i][j][2])
            {
            case 0: {
                xoc_split_schedule[i][j][0] = i * (new_k * w) + xoc_split_schedule[i][j][0];
                xoc_split_schedule[i][j][2] = 1;
                break;
            }
            case 1:
            case 4:
            case 5: {
                xoc_split_schedule[i][j][0] = i * (new_k * w) + xoc_split_schedule[i][j][0];
                break;
            }
            case 2: {
                xoc_split_schedule[i][j][2] = 3;
                break;
            }
            default:
                break;
            }
        }
    }

    // for(int i = 0; i < split_schedule.size(); i++) {
    //     print_smart_schedule(split_schedule[i]);
    //     printf("***************\n");
    // }
    // for(int i = 0; i < xoc_split_schedule.size(); i++) {
    //     print_comxor_schedule(xoc_split_schedule[i]);
    //     printf("-------------\n");
    // }

    int schedule_cpy = 0;
    int schedule_xor = 0;
    int common_cpy = 0;
    int common_xor = 0;
    for(int i = 0; i < split_schedule.size(); i++) {
        int j = 0;
        while(split_schedule[i][j][0] != -1) {
            if(split_schedule[i][j][4] == 0) {
                schedule_cpy++;
            }
            else {
                schedule_xor++;
            }
            j++;
        }
    }
    for(int i = 0; i < xoc_split_schedule.size(); i++) {
        for(int j = 0; j < xoc_split_schedule[i].size(); j++) {
            if(xoc_split_schedule[i][j][2] % 2 == 0) {
                common_cpy++;
            }
            else {
                common_xor++;
            }
        }
    }
    cout << "********************************************" << endl;
    cout << "split---schedlue copy nums:\t" << schedule_cpy << endl;
    cout << "split---schedlue xor nums:\t" << schedule_xor << endl;
    
    cout << "split---common copy nums:\t" << common_cpy << endl;
    cout << "split---common xor nums:\t" << common_xor << endl;
    


    // 执行调度操作

    // schedule_encode_single_chunk(char** data, int len, char**& parities, int K, int M, int W, int packetsize, int** en_schedule)
    // commonxor_encode_single_chunk(char **data, int len, char **&parities, int K, int M, int W, int packetsize, vector<int*> en_schedule, int intermedia_size)

    
    int loops = data_size / blocksize;
    char* dat = (char*)malloc(data_size);
    // cout << "***************************fenpeineicun1" << endl;
    cout << "**************************初始化数??" << endl;
    for(long long i = 0;i<data_size;i++)
    {
        dat[i] = rand();
    }
    char** par = malloc2d(m, data_size/k);
    char** parities = new char*[m];

    vector<vector<char*> > intermedias;
    for(int i = 0; i < intermedia_size.size(); i++) {
        vector<char*> intermedia;
        for(int j = 0;j<intermedia_size[i];j++)
        {
            char *tmp;
            posix_memalign((void**)&tmp,32,packetsize);
            intermedia.push_back(tmp);
        }
        intermedias.push_back(intermedia);
    }
    
    // cout << "***************************" << endl;

    // 公共异或吞吐量计??
    for(int i = 0;i<m;i++)
        parities[i] = par[i];
    // gettimeofday(&t0,NULL);
    for(int i = 0;i<loops;i++)
    {
        for(int j = 0; j < split_schedule.size(); j++) {
            commonxor_encode_single_chunk(dat+i*blocksize,blocksize,parities, k, m, w, packetsize, xoc_split_schedule[j], intermedias[j]);
        }
        for(int j = 0;j<m;j++)
        {
            parities[j] += blocksize/k;
        }
    }
    // gettimeofday(&t1,NULL);
    // diff = diff_us(t0,t1);
    printf(" !!! common xor------Encode cost %f us, size = %lld, speed = %f\n", T1 / 1000.0, data_size, (data_size / 1024.0/1024/1024 )/ (T1 /1000.0/1000 / 1000));
    // unsigned char* p = (unsigned char*)*par;
    // for(int i = 0; i < data_size / k * m; i++) {
    //     printf("%02x ", *p);
    //     p++;
    //     if((i + 1) % 20 == 0) {
    //         printf("\n");
    //     }
    // }
    // printf("\n");
    printf(" !!!xor_time: %f us, cpy_time: %f us\n", COM_XOR_TIME / 1000.0, COM_CPY_TIME / 1000.0);
    printf(" !!!all_time: %f us\n", (COM_XOR_TIME + COM_CPY_TIME) / 1000.0);
    printf("***************************************************\n");


    // 调度吞吐量计??
    for(int i = 0;i<m;i++)
        parities[i] = par[i];
    // gettimeofday(&t0,NULL);
    for(int i = 0;i<loops;i++)
    {
        for(int j = 0; j < split_schedule.size(); j++) {
            schedule_encode_single_chunk(dat+i*blocksize,blocksize,parities, k, m, w, packetsize, split_schedule[j]);
        }
        for(int j = 0;j<m;j++)
        {
            parities[j] += blocksize/k;
        }
    }
    // gettimeofday(&t1,NULL);
    // diff = diff_us(t0,t1);
    printf(" !!! schedule------Encode cost %f us, size = %lld, speed = %f\n", T0 / 1000.0, data_size,( data_size / 1024.0/1024/1024) / (T0 /1000.0/1000/1000));
    // unsigned char* p = (unsigned char*)*par;
    // for(int i = 0; i < data_size / k * m; i++) {
    //     printf("%02x ", *p);
    //     p++;
    //     if((i + 1) % 20 == 0) {
    //         printf("\n");
    //     }
    // }
    // printf("\n");
    long long xor_time = getXOR_TIME();
    long long cpy_time = getCPY_TIME();
    printf(" !!!xor_time: %f us, cpy_time: %f us\n", xor_time / 1000.0, cpy_time / 1000.0);
    printf(" !!!all_time: %f us\n", (xor_time + cpy_time) / 1000.0);
    
}

// 测试智能调度和公共异或，编码吞吐量的对比,10MB的数??
void test3(const int k, int m, const int w, const int nr , const int packetsize, long long data_size) {
    const int new_k = k / nr;
    const int blocksize = packetsize * k  *w;
    struct timeval t0,t1;
    long long  diff;
    while(data_size % blocksize != 0) {
        data_size++;
    }
    clock_t start, end;
    printf("**************k = %d, m = %d, w = %d, nr = %d, size = %d, packetsize = %d**************\n", k, m, w, nr, data_size, packetsize);

    int* cauchy = cauchy_good_general_coding_matrix(k, m, w);           // 生成编码矩阵
    jerasure_print_matrix(cauchy, m, k, w);

    int ones = get_sum_ones(cauchy, k, m , w);                           // 统计矩阵??1的个??
    cout << "bitmatrix ones: \t\t" << ones << endl;
    cout << "xor nums: \t\t\t" << ones - m * w << endl;

    int* bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, cauchy);     // 转换为位矩阵
    jerasure_print_bitmatrix(bitmatrix, m*w, k*w, w);
    
    vector<int*> split_bitmatrix = split(bitmatrix, k, m, w, nr);       // 将位矩阵进行分割
    
    
    // int xor_num = 0;                                                  // 统计使用智能调度后的异或总次??
    vector<int**> split_schedule;               // 保存智能调度，划分后的智能调??                     
    for(int i = 0; i < nr; i++) {
        int** schedule = jerasure_smart_bitmatrix_to_schedule(new_k, m, w, split_bitmatrix[i]);
        // xor_num += get_sum_schedules(schedule) - m * w;
        split_schedule.push_back(schedule);
    }
    // cout << "split---schedlue xor nums:\t" << xor_num + (nr - 1) * m * w << endl;


    // for(int i = 0; i < split_schedule.size(); i++) {
    //     print_smart_schedule(split_schedule[i]);
    //     printf("***************\n");
    // }

    // 恢复调度编号
    int ii = 0;
    while(split_schedule[0][ii][0] != -1) {
        split_schedule[0][ii][2] += (k - new_k);
        if(split_schedule[0][ii][0] >= new_k) {
            split_schedule[0][ii][0] += (k - new_k);
        }
        ii++;
    } 
    for(int i = 1; i < split_schedule.size(); i++) {
        int ii = 0;
        while(split_schedule[i][ii][0] != -1) {
            split_schedule[i][ii][2] += (k - new_k);
            split_schedule[i][ii][4] = 1;
            if(split_schedule[i][ii][0] >= new_k) {
                split_schedule[i][ii][0] += (k - new_k);
            }
            else {
                split_schedule[i][ii][0] = split_schedule[i][ii][0] + i * new_k; 
            }
            ii++;
        }
    }
    
    // for(int i = 0; i < split_schedule.size(); i++) {
    //     print_smart_schedule(split_schedule[i]);
    //     printf("***************\n");
    // }


    // for(int i = 0; i < split_schedule.size(); i++) {
    //     print_smart_schedule(split_schedule[i]);
    //     printf("***************\n");
    // }
    // for(int i = 0; i < xoc_split_schedule.size(); i++) {
    //     print_comxor_schedule(xoc_split_schedule[i]);
    //     printf("-------------\n");
    // }

    int schedule_cpy = 0;
    int schedule_xor = 0;
    int common_cpy = 0;
    int common_xor = 0;
    for(int i = 0; i < split_schedule.size(); i++) {
        int j = 0;
        while(split_schedule[i][j][0] != -1) {
            if(split_schedule[i][j][4] == 0) {
                schedule_cpy++;
            }
            else {
                schedule_xor++;
            }
            j++;
        }
    }
    cout << "********************************************" << endl;
    cout << "split---schedlue copy nums:\t" << schedule_cpy << endl;
    cout << "split---schedlue xor nums:\t" << schedule_xor << endl;

    // 执行调度操作

    // schedule_encode_single_chunk(char** data, int len, char**& parities, int K, int M, int W, int packetsize, int** en_schedule)
    // commonxor_encode_single_chunk(char **data, int len, char **&parities, int K, int M, int W, int packetsize, vector<int*> en_schedule, int intermedia_size)

    
    int loops = data_size / blocksize;
    char* dat = (char*)malloc(data_size);
    // cout << "***************************fenpeineicun1" << endl;
    cout << "**************************初始化数??" << endl;
    for(long long i = 0;i<data_size;i++)
    {
        dat[i] = rand();
    }
    char** par = malloc2d(m, data_size/k);
    char** parities = new char*[m];


    // 调度吞吐量计??
    for(int i = 0;i<m;i++)
        parities[i] = par[i];
    // gettimeofday(&t0,NULL);
    for(int i = 0;i<loops;i++)
    {
        for(int j = 0; j < split_schedule.size(); j++) {
            schedule_encode_single_chunk(dat+i*blocksize,blocksize,parities, k, m, w, packetsize, split_schedule[j]);
        }
        for(int j = 0;j<m;j++)
        {
            parities[j] += blocksize/k;
        }
    }
    // gettimeofday(&t1,NULL);
    // diff = diff_us(t0,t1);
    printf(" !!! schedule------Encode cost %f us, size = %lld, speed = %f\n", T0 / 1000.0, data_size,( data_size / 1024.0/1024/1024) / (T0 /1000.0/1000/1000));
    // unsigned char* p = (unsigned char*)*par;
    // for(int i = 0; i < data_size / k * m; i++) {
    //     printf("%02x ", *p);
    //     p++;
    //     if((i + 1) % 20 == 0) {
    //         printf("\n");
    //     }
    // }
    // printf("\n");
    long long xor_time = getXOR_TIME();
    long long cpy_time = getCPY_TIME();
    printf(" !!!xor_time: %f us, cpy_time: %f us\n", xor_time / 1000.0, cpy_time / 1000.0);
    printf(" !!!all_time: %f us\n", (xor_time + cpy_time) / 1000.0);
    
}

//**************************************************************************
//    同时发现多个异或的测??
const int TEST4_K = 10;
const int TEST4_M = 4;
const int TEST4_W = 8;
const int BIT_SIZE = TEST4_K * TEST4_W;
// const int SELECT_N_AND = 2;
template<size_t N>
struct bitsetAndSet {
    set<int> sets;
    size_t oneNums;
    bitset<N> andSet;
};
template<size_t N>
bool operator<(const bitsetAndSet<N>& b1, const bitsetAndSet<N>& b2) {
    return b1.oneNums > b2.oneNums;
}

template<size_t N>
vector<bitset<N>> bitmatrix_to_bitsets(int* bitmatrix, int k, int m, int w) {
    vector<bitset<N>> bitsets;
    for(int i = 0; i < m * w; i++) {
        int start = i * k * w;
        string s;
        for(int j = 0; j < k * w; j++) {
            s.push_back(bitmatrix[start + j] + '0');
            // s += (bitmatrix[start + j] + "");
        }
        bitset<N> bs(s);
        bitsets.push_back(bs);
    }
    return bitsets;
}

ostream& operator<<(ostream& os, set<int>& s) {
    for(auto a: s) {
        cout << a << " ";
    }
    return os;
}

template<size_t N>
ostream& operator<<(ostream& os, vector<bitsetAndSet<N>>& vc) {
    for(auto st: vc) {
        cout << st.sets << "****" << st.oneNums << "****" << st.andSet << endl;
    }
    return os;
}
template<class T>
ostream& operator<<(ostream& os, vector<T>& vc) {
    for(auto c: vc) {
        cout << c << endl;
    }
    return os;
}

// 统计所有n个packet相与的结??
template<size_t N>
void operator_and_N(vector<bitset<N>>& bitsets, int start, const int n, 
vector<bitsetAndSet<N>>& result, set<int>& set_path, vector<bitset<N>>& and_path) {
    // if(start + n > bitsets.size()) {
    //     return;
    // }
    if(set_path.size() == n) {
        bitsetAndSet<N> path {set_path, and_path.back().count(), and_path.back()};
        result.push_back(path);
        // and_path.pop_back();
        return;
    }
    
    for(int i = start; i < bitsets.size(); i++) {
        set_path.insert(i);
        and_path.push_back(and_path.back());
        // cout << "*********************" << and_path.back() << endl;
        and_path.back() &= bitsets[i];
        // cout << "*********************" << bitsets[i] << "****" << i << endl;
        // cout << set_path << "*******************" << endl;
        // for(auto a: and_path) {
        //     cout << a << "*******************" << endl;
        // }
        operator_and_N(bitsets, i + 1, n, result, set_path, and_path);
        set_path.erase(i);
        and_path.pop_back();
        // cout << set_path << "*******************" << endl;
        // for(auto a: and_path) {
        //     cout << a << "*******************" << endl;
        // }
    }
}

// 求出“最优”的中间结果，并更改原先的位矩阵
// 返回因为公共异或而减少的异或次数
template<size_t N>
int get_optimation_update_bitsets(vector<bitsetAndSet<N>> result, vector<bitset<N>>& bitsets, vector<bitsetAndSet<N>>& intermedia_result) {
    set<int> hash_set;
    auto it = result.begin();
    int desc = 0;
    int and_size = result[0].sets.size();           // 与运算的packet的数??
    while(hash_set.size() != (bitsets.size())) {
        // cout << it->oneNums  << endl;
        if(it->oneNums < 2) {
            break;
        }
        auto setIt = it->sets.begin();
        while(setIt != it->sets.end()) {
            if(hash_set.find(*setIt) != hash_set.end()) {
                break;
            }
            ++setIt;
        }
        if(setIt != it->sets.end()) {
            ++it;
            continue;
        }
        else {
            desc += (it->oneNums - 1) * (and_size - 1);
            for(auto a: it->sets) {
                hash_set.insert(a);             // 将元素插入集合，代表已经选取
                bitsets[a] ^= it->andSet;       // 将已经被选取的packetId置为0
            }
            intermedia_result.push_back(*it);
            ++it;
        }
    }
    return desc;
}
// 执行一次算法，统计减小的异或次??, n_and为参与与运算的packet数目
template<size_t N>
int do_operate(vector<bitset<N>>& bitsets, int n_and, int& inc_copy) {
    vector<bitsetAndSet<BIT_SIZE>> result;
    set<int> set_path;
    vector<bitset<BIT_SIZE>> and_path;
    bitset<BIT_SIZE> bs;
    bs.set();
    and_path.push_back(bs);
    vector<bitsetAndSet<BIT_SIZE>> intermedia_result;
    
    operator_and_N<BIT_SIZE>(bitsets, 0, n_and,  result, set_path, and_path);
    sort(result.begin(), result.end());
    // for(auto r: result) {
    //     cout << r.sets << "****" << r.oneNums << "****" << r.andSet << endl;
    // }
    // cout << result << endl << "*******************" << endl;
    
    int desc = get_optimation_update_bitsets(result, bitsets, intermedia_result);
    cout << intermedia_result << endl;
    inc_copy = intermedia_result.size();
    return desc;
}


void test_cauchy(const int k, const int m, const int w) {
    printf("*****************************************k = %d, m = %d, w = %d*****************************************\n", k, m, w);

    int* cauchy = cauchy_original_coding_matrix(k, m, w);       // 生成编码矩阵
    // int* cauchy = cauchy_good_general_coding_matrix(k, m, w);       // 生成编码矩阵
    
    jerasure_print_matrix(cauchy, m, k, w);
    int ones = get_sum_ones(cauchy, k, m , w);                           // 统计矩阵??1的个??
    int xors = ones - w * m;                                             // 异或次数


    int* bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, cauchy);     // 转换为位矩阵
    // jerasure_print_bitmatrix(bitmatrix, m*w, k*w, w);

    int** schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
    // print_smart_schedule(schedule);
    int schedule_copy = get_copy_schedules(schedule);
    int schedule_xor = get_xor_schedules(schedule);
    

    ZOXC oxc(k*w, m*w);
    oxc.grouping_1s(bitmatrix, false);
    int common_copy = get_copy_common(oxc.schedule);
    int common_xor = get_xor_common(oxc.schedule);
    cout << oxc.intermedia_schedule.size() << endl;
    // print_comxor_schedule(oxc.intermedia_schedule);

    auto bitsets = bitmatrix_to_bitsets<BIT_SIZE>(bitmatrix, k, m, w);


    int mine_copy = m * w;
    int mine_xor = xors;
    int desc = 0;
    int inc_copy;
    vector<int> lists{2,2,2,2,4};
    for(int i = 0; i < lists.size(); i++) {
        desc = do_operate(bitsets, lists[i], inc_copy);
        mine_copy += inc_copy;
        mine_xor -= desc;
        cout << mine_copy << "--------------------" << mine_xor << endl;
    }
    cout << "bitmatrix ones: \t\t" << ones << endl;
    cout << "xor nums: \t\t\t" << ones - m * w << endl;

    cout << "********************************************" << endl;
    cout << "schedule copy nums:\t" << schedule_copy << endl;
    cout << "common copy nums:\t" << common_copy << endl;
    cout << "mine copy nums:\t" << mine_copy << endl << endl;

    cout << "schedule xor nums:\t" << schedule_xor << endl;
    cout << "common xor nums:\t" << common_xor << endl;
    cout << "mine xor nums:\t" << mine_xor  << endl << endl;

    cout << "schedule sum nums:\t" << schedule_copy + schedule_xor << endl;
    cout << "common sum nums:\t" << common_copy + common_xor << endl;
    cout << "mine sum nums:\t" << mine_copy + mine_xor << endl << endl;
    printf("*****************************************k = %d, m = %d, w = %d", k, m, w);
    cout << "***********";
    for(auto a: lists) {
        cout << a << " ";
    }
    cout << "*******************" << endl << endl;
    
}


void test_cauchy_good(const int k, const int m, const int w) {
    printf("*****************************************k = %d, m = %d, w = %d*****************************************\n", k, m, w);

    int* cauchy = cauchy_good_general_coding_matrix(k, m, w);       // 生成编码矩阵
    // int* cauchy = cauchy_good_general_coding_matrix(k, m, w);       // 生成编码矩阵
    
    // jerasure_print_matrix(cauchy, m, k, w);
    int ones = get_sum_ones(cauchy, k, m , w);                           // 统计矩阵??1的个??
    int xors = ones - w * m;                                             // 异或次数


    int* bitmatrix = jerasure_matrix_to_bitmatrix(k, m, w, cauchy);     // 转换为位矩阵
    // jerasure_print_bitmatrix(bitmatrix, m*w, k*w, w);

    int** schedule = jerasure_smart_bitmatrix_to_schedule(k, m, w, bitmatrix);
    int schedule_copy = get_copy_schedules(schedule);
    int schedule_xor = get_xor_schedules(schedule);
    

    ZOXC oxc(k*w, m*w);
    oxc.grouping_1s(bitmatrix, false);
    int common_copy = get_copy_common(oxc.schedule);
    int common_xor = get_xor_common(oxc.schedule);
    cout << oxc.intermedia_schedule.size() << endl;

    int* new_cauchy = cauchy + k;
    int new_m = m - 1;
    // jerasure_print_matrix(new_cauchy, new_m, k, w);
    int* new_bitmatrix = jerasure_matrix_to_bitmatrix(k, new_m, w, new_cauchy); 

    auto bitsets = bitmatrix_to_bitsets<BIT_SIZE>(new_bitmatrix, k, new_m, w);
    // for(auto bits: bitsets) {
    //     cout << bits << endl;
    // }
   


    int mine_copy = m * w;
    int mine_xor = xors;
    int desc = 0;
    int inc_copy;
    vector<int> lists{4,4,3,3,2,2,2};
    for(int i = 0; i < lists.size(); i++) {
        desc = do_operate(bitsets, lists[i], inc_copy);
        mine_copy += inc_copy;
        mine_xor -= desc;
        cout << mine_copy << "--------------------" << mine_xor << endl;
    }
    cout << "bitmatrix ones: \t\t" << ones << endl;
    cout << "xor nums: \t\t\t" << ones - m * w << endl;

    cout << "********************************************" << endl;
    cout << "schedule copy nums:\t" << schedule_copy << endl;
    cout << "common copy nums:\t" << common_copy << endl;
    cout << "mine copy nums:\t" << mine_copy << endl << endl;

    cout << "schedule xor nums:\t" << schedule_xor << endl;
    cout << "common xor nums:\t" << common_xor << endl;
    cout << "mine xor nums:\t" << mine_xor  << endl << endl;

    cout << "schedule sum nums:\t" << schedule_copy + schedule_xor << endl;
    cout << "common sum nums:\t" << common_copy + common_xor << endl;
    cout << "mine sum nums:\t" << mine_copy + mine_xor << endl << endl;
    printf("*****************************************k = %d, m = %d, w = %d", k, m, w);
    cout << "***********";
    for(auto a: lists) {
        cout << a << " ";
    }
    cout << "*******************" << endl << endl;
    
}


//**************************************************************************



int main(int argc, char* argv[]) {
    // const int k = atoi(argv[1]);
    // const int m = atoi(argv[2]);
    // const int w = atoi(argv[3]);
    test_cauchy(TEST4_K, TEST4_M, TEST4_W);
    return 0;
}




/*
main.cpp
Tianli Zhou

Fast Erasure Coding for Data Storage: A Comprehensive Study of the Acceleration Techniques

Revision 1.0
Mar 20, 2019

Tianli Zhou
Department of Electrical & Computer Engineering
Texas A&M University
College Station, TX, 77843
zhoutianli01@tamu.edu

Copyright (c) 2019, Tianli Zhou
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in
  the documentation and/or other materials provided with the
  distribution.

- Neither the name of the Texas A&M University nor the names of its
  contributors may be used to endorse or promote products derived
  from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/
/*

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include "Example/zexample.h"

void usage()
{
    printf("Usage:\n");
    printf("\t./zerasure test_cost_weight [size] [loops]\n");
    printf("\t\t- Test the weight of memcpy and XOR in the schedule, \n\t\tusing data size [size] and run [loops] times\n\n");

    printf("\t... | ./zerasure single cost_func[0..2] strategy[0..7]\n");
    printf("\t\t- Read K,M,W,X,Y from stdin, generate schedule for given X,Y array\n");
    printf("\t\tusing cost function [cost_func] and strategy [strategy]\n\n");

    printf("\t./zerasure sa K M W S acc_rate cost_func[0..2] strategy[0..7]\n");
    printf("\t\t- Perform simulated annealing optimization, stop if minimum\n");
    printf("\t\thaven't change for S consecutive iterations. Initial rate\n");
    printf("\t\tof accept worse successor is given as acc_rate\n");
    printf("\t\tusing cost function [cost_func] and strategy [strategy]\n\n");

    printf("\t./zerasure ge K M W S init_population select_rate crossover_rate\n");
    printf("\t           tmutation_rate max_population cost_func[0..2] strategy[0..7]\n");
    printf("\t\t- Perform genetic optimization using given parameters, stop if minimum\n");
    printf("\t\thaven't change for S consecutive iterations. Initial rate\n");
    printf("\t\tusing cost function [cost_func] and strategy [strategy]\n\n");

    printf("\t... | ./zerasure code packetsize strategy[0..7]\n");
    printf("\t\t- Read K,M,W,X,Y from stdin, generate schedule for given\n");
    printf("\t\tusing cost function [cost_func] and strategy [strategy]\n");
}

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        usage();
        exit(-1);
    }

#ifdef VEC128
    printf("Using 128-bit vectorization\n");
#elif VEC256
    printf("Using 256-bit vectorization\n");
#else
    printf("Unknown Vectorization\n");
#endif
    printf("\n\n ********* Start zerasure *********\n\n");

    if(strcmp(argv[1], "test_cost_weight") == 0)
        ZExample::test_cost_weight(argc,argv);

    if(strcmp(argv[1], "single") == 0)
        ZExample::single(argc,argv);

    if(strcmp(argv[1], "sa") == 0)
        ZExample::sa(argc,argv);

    if(strcmp(argv[1], "ge") == 0)
        ZExample::ge(argc,argv);

    if(strcmp(argv[1], "code") == 0)
        ZExample::code(argc,argv);

    printf("\n\n *********  End zerasure **********\n\n");
    return 0;
}
*/