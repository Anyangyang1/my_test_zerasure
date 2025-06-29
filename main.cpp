
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
#include<cmath>
using namespace std;
const int K = 48;
const int M = 4;
const int W = 8;
const int FIRST_SELECT = 6;

void print_ids(int* dm_ids) {
	printf("************ids*************");
	int i = 0;
	for(i = 0; i < K; i++) {
		printf("%d ", dm_ids[i]);
	}
	printf("\n");
}
// 递归函数，用于生成组合
void generate_combinations(int start, int n, int k, int* arr, int* current, int current_size, int** result, int* count) {
    if (current_size == k) {
        // 保存当前组合到结果数组中
        result[*count] = (int*)malloc(k * sizeof(int));
        for (int i = 0; i < k; i++) {
            result[*count][i] = current[i];
        }
        (*count)++;
        return;
    }

    for (int i = start; i < n; i++) {
        current[current_size] = arr[i];
        generate_combinations(i + 1, n, k, arr, current, current_size + 1, result, count);
    }
}

// 主函数，用于生成所有 k 组合
int** get_combinations(int n, int k, int* arr, int* count) {
    // 计算组合的总数
    int total_combinations = 1;
    for (int i = 0; i < k; i++) {
        total_combinations *= (n - i);
        total_combinations /= (i + 1);
    }

    // 分配内存来存储所有组合
    int** result = (int**)malloc(total_combinations * sizeof(int*));
    int* current = (int*)malloc(k * sizeof(int));
    *count = 0;

    // 生成组合
    generate_combinations(0, n, k, arr, current, 0, result, count);

    // 释放临时数组
    free(current);

    return result;
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

// 递归函数，用于生成组合
void generate_combinations2(int start, int n, int k, std::vector<int>& current, std::vector<std::vector<int>>& result, set<int>& isExist) {
    if (current.size() == k) {
        result.push_back(current);
        return;
    }

    for (int i = start; i <= n; ++i) {
        if(isExist.count(i % k) > 0) {
            continue;
        }
        current.push_back(i);
        isExist.insert(i % k);
        generate_combinations2(i + 1, n, k, current, result, isExist);
        current.pop_back();
        isExist.erase(i % k);
    }
}

// 主函数，用于生成所有 k 组合
std::vector<std::vector<int>> get_combinations2(int n, int k) {
    std::vector<std::vector<int>> result;
    std::vector<int> current;
    set<int> isExist;
    generate_combinations2(0, n - 1, k, current, result, isExist);
    return result;
}

template<class T>
ostream& operator<<(ostream& os, vector<T>& vc) {
    for(auto c: vc) {
        cout << c << " ";
    }
    return os;
}


// 高斯消元法计算 01 矩阵的秩
int calculate_matrix_rank(std::vector<std::vector<int>>& matrix) {
    int n = matrix.size();  // 矩阵的行数
    int m = matrix[0].size();  // 矩阵的列数
    int rank = 0;  // 矩阵的秩

    for (int col = 0; col < m; ++col) {
        // 查找当前列中第一个为 1 的行
        int pivot_row = -1;
        for (int row = rank; row < n; ++row) {
            if (matrix[row][col] == 1) {
                pivot_row = row;
                break;
            }
        }

        // 如果找到这样的行，进行行变换
        if (pivot_row != -1) {
            // 交换当前行和找到的行
            std::swap(matrix[rank], matrix[pivot_row]);

            // 使用当前行消去其他行的第 col 列
            for (int row = 0; row < n; ++row) {
                if (row != rank && matrix[row][col] == 1) {
                    for (int j = 0; j < m; ++j) {
                        matrix[row][j] ^= matrix[rank][j];
                    }
                }
            }

            // 增加秩
            ++rank;
        }
    }

    return rank;
}



void force_research() {
    // 生成编码矩阵
    int* matrix = cauchy_original_coding_matrix(K, M, W);
	jerasure_print_matrix(matrix, M, K, W);
	int* bitmatrix = jerasure_matrix_to_bitmatrix(K, M, W, matrix);
	jerasure_print_bitmatrix(bitmatrix, M * W, K * W, W);

	// 生成所有组合的修复方式
	int i, j;
    int count = 0;
    int* nums = (int*) malloc (sizeof(int) * (K + M - 1));
    for(i = 1; i <= K + M - 1; i++) {
        nums[i - 1] = i;
    }
    // 分配内存来存储所有组合
    int **combinations = get_combinations(K + M - 1, K, nums, &count);
    for (i = 0; i < count; i++) {
        for (j = 0; j < K; j++) {
            printf("%d ", combinations[i][j]);
        }
        printf("\n");
        // free(combinations[i]);  // 释放每个组合的内存
    }
    // free(combinations);  // 释放组合数组的内存
    // free(nums);
    

    int decode_matrix[K * K];
	int erased[K + M] ;
	int dm_ids[K];
    int all_decode[count * (K + M)];
    for(i = 0; i < count; i++) {
        for(j = 0; j < K + M; j++) {
            all_decode[i * (K + M) + j] = 0;
        }
    }
    for(i = 0; i < count; i++) {
		for(j = 0; j < K + M; j++) {
			erased[j] = 1;
		}
		for(j = 0; j < K; j++) {
			erased[combinations[i][j]] = 0;
		}
		jerasure_make_decoding_matrix(K, M, W, matrix, erased, decode_matrix, dm_ids);
        print_ids(dm_ids);
		jerasure_print_matrix(decode_matrix, K, K, W);
		for(j = 0; j < K; j++) {
            all_decode[i * (K + M) + dm_ids[j]] = decode_matrix[j];
        }
		
	}
    for(i = 0; i < count; i++) {
        for(j = 0; j < K + M; j++) {
            printf("%d\t", all_decode[i * (K + M) + j]);
        }
        printf("\n");
    }
    int* all_bitmatrix = jerasure_matrix_to_bitmatrix(K + M, count ,W, all_decode);
    jerasure_print_bitmatrix(all_bitmatrix, count*W, (K+M)*W, W);
    cout << " *************************" << endl;

    int minV = INT32_MAX;
    vector<vector<int>> gK = get_combinations2(count*W, W);
    for(vector<int> g: gK) {
        cout << g << "***********";
        int sum = 0;
        for(int col = 1; col < K + M; col++) {
            vector<vector<int>> nums1;
            for(int row: g) {
                vector<int> nums2;
                for(int k = 0; k < W; k++) {
                    nums2.push_back(all_bitmatrix[row * (K + M) * W + col * W + k]);
                }
                nums1.push_back(nums2);
            }
            int rank = calculate_matrix_rank(nums1);
            cout << rank << " ";
            sum += rank;
        }
        cout << "************";
        cout << sum;
        minV = min(minV, sum);
        cout << "**********" << minV << endl;
    }
    cout << "minV: " << minV << endl;
}

int getNum(int* bit, int w) {
    int ans = 0;
    for(int i = 0; i < w; i++) {
        ans = (ans << 1) + bit[i];
    }
    return ans;
}

void transform_bit_to_integer(int* all_bitmatrix, int row, int col, int w, vector<vector<int>>& all_recovery) {
    for(int i = 0; i < row; i++) {
        vector<int> tmp;
        for(int j = 0; j < col; j++) {
            tmp.push_back(getNum(all_bitmatrix + (i * col * w) + j * w, w));
        }
        all_recovery.push_back(tmp);
    }
}

int getSimilar(vector<int>& nums, vector<set<int>> xor_hash) {
    int ans = 0;
    for(int i = 0; i < nums.size(); i++) {
        if(nums[i] != 0 && xor_hash[i].count(nums[i]) > 0) {
            ans++;
        }
    }
    return ans;
}


int getRowNum(vector<vector<int>> &all_recovery, vector<set<int>> xor_hash, set<int> exist_hash, int w, int &ans) {
    int maxV = -1;
    for(int i = 0; i < all_recovery.size(); i++) {
        if(exist_hash.count(i % w) == 0) {
            int sim = getSimilar(all_recovery[i], xor_hash);
            if(sim > maxV) {
                maxV = sim;
                ans = i;
            }
        }
    }
    return maxV;
}


void add_xor_hash(set<int>& hash, int val) {
    vector<int> tmp;
    for(int h: hash) {
        tmp.push_back(h);
    }
    for(int t: tmp) {
        hash.insert(t ^ val);
    }
    hash.insert(val);
}


void tan_heart() {
    int* matrix = cauchy_original_coding_matrix(K, M, W);
	// jerasure_print_matrix(matrix, M, K, W);
	int* bitmatrix = jerasure_matrix_to_bitmatrix(K, M, W, matrix);
	// jerasure_print_bitmatrix(bitmatrix, M * W, K * W, W);

	// 生成所有组合的修复方式
	int i, j;
    int count = 0;
    int* nums = (int*) malloc (sizeof(int) * (K + M - 1));
    for(i = 1; i <= K + M - 1; i++) {
        nums[i - 1] = i;
    }
    // 分配内存来存储所有组合
    int **combinations = get_combinations(K + M - 1, K, nums, &count);
    // for (i = 0; i < count; i++) {
    //     for (j = 0; j < K; j++) {
    //         printf("%d ", combinations[i][j]);
    //     }
    //     printf("\n");
    //     // free(combinations[i]);  // 释放每个组合的内存
    // }
    // free(combinations);  // 释放组合数组的内存
    // free(nums);
    

    int decode_matrix[K * K];
	int erased[K + M] ;
	int dm_ids[K];
    int all_decode[count * (K + M)];
    for(i = 0; i < count; i++) {
        for(j = 0; j < K + M; j++) {
            all_decode[i * (K + M) + j] = 0;
        }
    }
    for(i = 0; i < count; i++) {
		for(j = 0; j < K + M; j++) {
			erased[j] = 1;
		}
		for(j = 0; j < K; j++) {
			erased[combinations[i][j]] = 0;
		}
		jerasure_make_decoding_matrix(K, M, W, matrix, erased, decode_matrix, dm_ids);
        // print_ids(dm_ids);
		// jerasure_print_matrix(decode_matrix, K, K, W);
		for(j = 0; j < K; j++) {
            all_decode[i * (K + M) + dm_ids[j]] = decode_matrix[j];
        }
		
	}
    // for(i = 0; i < count; i++) {
    //     for(j = 0; j < K + M; j++) {
    //         printf("%d\t", all_decode[i * (K + M) + j]);
    //     }
    //     printf("\n");
    // }
    int* all_bitmatrix = jerasure_matrix_to_bitmatrix(K + M, count ,W, all_decode);
    // cout << " *************************" << endl;
    // jerasure_print_bitmatrix(all_bitmatrix, count*W, (K+M)*W, W);
    // cout << " *************************" << endl;

    vector<vector<int>> all_recovery;
    transform_bit_to_integer(all_bitmatrix, count*W, K + M, W, all_recovery);
    // for(vector<int> tmp: all_recovery) {
    //     cout << tmp << endl;
    // }

    vector<set<int>> xor_hash(K + M, set<int>());
    set<int> exist_hash;
    
    for(i = 0; i < K + M; i++) {
        if(all_recovery[FIRST_SELECT][i] != 0) {
            xor_hash[i].insert(all_recovery[FIRST_SELECT][i]);
        }
    }
    exist_hash.insert(FIRST_SELECT % W);
    int canDecrease = 0;
    vector<int> collect_idx;
    collect_idx.push_back(FIRST_SELECT);
    cout << "decrease: ";
    while(exist_hash.size() != W) {
        int idx;
        int decrease = getRowNum(all_recovery, xor_hash, exist_hash, W, idx);
        cout << decrease << " ";
        canDecrease += decrease;
        for(i = 0; i < K + M; i++) {
            if(all_recovery[idx][i] != 0) {
                add_xor_hash(xor_hash[i], all_recovery[idx][i]);
            }
        }
        exist_hash.insert(idx % W);
        collect_idx.push_back(idx);
    }
    cout << "**************K, M, W*************" << K << "," << M << "," << W << endl;
    cout << "origin_blocks: " << K * W << endl;
    cout << "my-tra:" << K * W - canDecrease << endl;

    cout << "collect_idx: " << collect_idx << endl;
    cout << "first_select: " << FIRST_SELECT << endl;
    // for(auto xh: xor_hash) {
    //     for(auto x: xh) {
    //         cout << x << " ";
    //     }
    //     cout << endl;
    // }

}


void tan_heart_combine_reserve() {
    int* matrix = cauchy_original_coding_matrix(K, M, W);
	// jerasure_print_matrix(matrix, M, K, W);
	int* bitmatrix = jerasure_matrix_to_bitmatrix(K, M, W, matrix);
	// jerasure_print_bitmatrix(bitmatrix, M * W, K * W, W);

	// 生成所有组合的修复方式
	int i, j;
    int count = 0;
    int* nums = (int*) malloc (sizeof(int) * (K + M - 1));
    for(i = 1; i <= K + M - 1; i++) {
        nums[i - 1] = i;
    }
    // 分配内存来存储所有组合
    int **combinations = get_combinations(K + M - 1, M - 1, nums, &count);
    // for (i = 0; i < count; i++) {
    //     for (j = 0; j < M - 1; j++) {
    //         printf("%d ", combinations[i][j]);
    //     }
    //     printf("\n");
    //     // free(combinations[i]);  // 释放每个组合的内存
    // }
    // // free(combinations);  // 释放组合数组的内存
    // // free(nums);
    

    int decode_matrix[K * K];
	int erased[K + M] ;
	int dm_ids[K];
    int all_decode[count * (K + M)];
    for(i = 0; i < count; i++) {
        for(j = 0; j < K + M; j++) {
            all_decode[i * (K + M) + j] = 0;
        }
    }
    for(i = 0; i < count; i++) {
		for(j = 0; j < K + M; j++) {
			erased[j] = 0;
		}
		for(j = 0; j < M - 1; j++) {
			erased[combinations[i][j]] = 1;
		}
        erased[0] = 1;
		jerasure_make_decoding_matrix(K, M, W, matrix, erased, decode_matrix, dm_ids);
        // print_ids(dm_ids);
		// jerasure_print_matrix(decode_matrix, K, K, W);
		for(j = 0; j < K; j++) {
            all_decode[i * (K + M) + dm_ids[j]] = decode_matrix[j];
        }
		
	}
    // for(i = 0; i < count; i++) {
    //     for(j = 0; j < K + M; j++) {
    //         printf("%d\t", all_decode[i * (K + M) + j]);
    //     }
    //     printf("\n");
    // }
    int* all_bitmatrix = jerasure_matrix_to_bitmatrix(K + M, count ,W, all_decode);
    // cout << " *************************" << endl;
    // jerasure_print_bitmatrix(all_bitmatrix, count*W, (K+M)*W, W);
    // cout << " *************************" << endl;

    vector<vector<int>> all_recovery;
    transform_bit_to_integer(all_bitmatrix, count*W, K + M, W, all_recovery);
    // for(vector<int> tmp: all_recovery) {
    //     cout << tmp << endl;
    // }

    vector<set<int>> xor_hash(K + M, set<int>());
    set<int> exist_hash;
    
    for(i = 0; i < K + M; i++) {
        if(all_recovery[FIRST_SELECT][i] != 0) {
            xor_hash[i].insert(all_recovery[FIRST_SELECT][i]);
        }
    }
    exist_hash.insert(FIRST_SELECT % W);
    int canDecrease = 0;
    vector<int> collect_idx;
    collect_idx.push_back(FIRST_SELECT);
    cout << "decrease: ";
    while(exist_hash.size() != W) {
        int idx;
        int decrease = getRowNum(all_recovery, xor_hash, exist_hash, W, idx);
        cout << decrease << " ";
        canDecrease += decrease;
        for(i = 0; i < K + M; i++) {
            if(all_recovery[idx][i] != 0) {
                add_xor_hash(xor_hash[i], all_recovery[idx][i]);
            }
        }
        exist_hash.insert(idx % W);
        collect_idx.push_back(idx);
    }
    cout << "**************K, M, W*************" << K << "," << M << "," << W << endl;
    cout << "origin_blocks: " << K * W << endl;
    cout << "my-tra:" << K * W - canDecrease << endl;

    cout << "collect_idx: " << collect_idx << endl;
    cout << "first_select: " << FIRST_SELECT << endl;
    // for(auto xh: xor_hash) {
    //     for(auto x: xh) {
    //         cout << x << " ";
    //     }
    //     cout << endl;
    // }
    int sum = 0;
    cout << "rank: ";
    for(int col = 1; col < K + M; col++) {
        vector<vector<int>> nums1;
        for(int row: collect_idx) {
            vector<int> nums2;
            for(int k = 0; k < W; k++) {
                nums2.push_back(all_bitmatrix[row * (K + M) * W + col * W + k]);
            }
            nums1.push_back(nums2);
        }
        int rank = calculate_matrix_rank(nums1);
        cout << rank << " ";
        sum += rank;
    }
    cout << endl;
    cout << "sum: " << sum << endl;
}


void tan_heart_with_print_combine_reserve() {
    int* matrix = cauchy_original_coding_matrix(K, M, W);
	jerasure_print_matrix(matrix, M, K, W);
	int* bitmatrix = jerasure_matrix_to_bitmatrix(K, M, W, matrix);
	jerasure_print_bitmatrix(bitmatrix, M * W, K * W, W);

	// 生成所有组合的修复方式
	int i, j;
    int count = 0;
    int* nums = (int*) malloc (sizeof(int) * (K + M - 1));
    for(i = 1; i <= K + M - 1; i++) {
        nums[i - 1] = i;
    }
    // 分配内存来存储所有组合
    int **combinations = get_combinations(K + M - 1, M - 1, nums, &count);
    free(nums);
    for (i = 0; i < count; i++) {
        for (j = 0; j < M - 1; j++) {
            printf("%d ", combinations[i][j]);
        }
        printf("\n");
        // free(combinations[i]);  // 释放每个组合的内存
    }
    // free(combinations);  // 释放组合数组的内存
    // free(nums);
    

    int decode_matrix[K * K];
	int erased[K + M] ;
	int dm_ids[K];
    int all_decode[count * (K + M)];
    for(i = 0; i < count; i++) {
        for(j = 0; j < K + M; j++) {
            all_decode[i * (K + M) + j] = 0;
        }
    }
    for(i = 0; i < count; i++) {
		for(j = 0; j < K + M; j++) {
			erased[j] = 0;
		}
		for(j = 0; j < M - 1; j++) {
			erased[combinations[i][j]] = 1;
		}
        free(combinations[i]); 
        erased[0] = 1;
		jerasure_make_decoding_matrix(K, M, W, matrix, erased, decode_matrix, dm_ids);
        print_ids(dm_ids);
		jerasure_print_matrix(decode_matrix, K, K, W);
		for(j = 0; j < K; j++) {
            all_decode[i * (K + M) + dm_ids[j]] = decode_matrix[j];
        }
	}
    free(combinations);
    for(i = 0; i < count; i++) {
        for(j = 0; j < K + M; j++) {
            printf("%d\t", all_decode[i * (K + M) + j]);
        }
        printf("\n");
    }
    int* all_bitmatrix = jerasure_matrix_to_bitmatrix(K + M, count ,W, all_decode);
    cout << " *************************" << endl;
    jerasure_print_bitmatrix(all_bitmatrix, count*W, (K+M)*W, W);
    cout << " *************************" << endl;

    vector<vector<int>> all_recovery;
    transform_bit_to_integer(all_bitmatrix, count*W, K + M, W, all_recovery);
    for(vector<int> tmp: all_recovery) {
        cout << tmp << endl;
    }

    vector<set<int>> xor_hash(K + M, set<int>());
    set<int> exist_hash;
    
    for(i = 0; i < K + M; i++) {
        if(all_recovery[FIRST_SELECT][i] != 0) {
            xor_hash[i].insert(all_recovery[FIRST_SELECT][i]);
        }
    }
    exist_hash.insert(FIRST_SELECT % W);
    int canDecrease = 0;
    vector<int> collect_idx;
    collect_idx.push_back(FIRST_SELECT);
    cout << "decrease: ";
    while(exist_hash.size() != W) {
        int idx;
        int decrease = getRowNum(all_recovery, xor_hash, exist_hash, W, idx);
        cout << decrease << " ";
        canDecrease += decrease;
        for(i = 0; i < K + M; i++) {
            if(all_recovery[idx][i] != 0) {
                add_xor_hash(xor_hash[i], all_recovery[idx][i]);
            }
        }
        exist_hash.insert(idx % W);
        collect_idx.push_back(idx);
    }
    cout << "**************K, M, W*************" << K << "," << M << "," << W << endl;
    cout << "origin_blocks: " << K * W << endl;
    cout << "my-tra:" << K * W - canDecrease << endl;

    cout << "collect_idx: " << collect_idx << endl;
    cout << "first_select: " << FIRST_SELECT << endl;
    for(auto xh: xor_hash) {
        for(auto x: xh) {
            cout << x << " ";
        }
        cout << endl;
    }

}



void tan_heart_with_print() {
    int* matrix = cauchy_original_coding_matrix(K, M, W);
	jerasure_print_matrix(matrix, M, K, W);
	int* bitmatrix = jerasure_matrix_to_bitmatrix(K, M, W, matrix);
	jerasure_print_bitmatrix(bitmatrix, M * W, K * W, W);

	// 生成所有组合的修复方式
	int i, j;
    int count = 0;
    int* nums = (int*) malloc (sizeof(int) * (K + M - 1));
    for(i = 1; i <= K + M - 1; i++) {
        nums[i - 1] = i;
    }
    // 分配内存来存储所有组合
    int **combinations = get_combinations(K + M - 1, K, nums, &count);
    for (i = 0; i < count; i++) {
        for (j = 0; j < K; j++) {
            printf("%d ", combinations[i][j]);
        }
        printf("\n");
        // free(combinations[i]);  // 释放每个组合的内存
    }
    // free(combinations);  // 释放组合数组的内存
    // free(nums);
    

    int decode_matrix[K * K];
	int erased[K + M] ;
	int dm_ids[K];
    int all_decode[count * (K + M)];
    for(i = 0; i < count; i++) {
        for(j = 0; j < K + M; j++) {
            all_decode[i * (K + M) + j] = 0;
        }
    }
    for(i = 0; i < count; i++) {
		for(j = 0; j < K + M; j++) {
			erased[j] = 1;
		}
		for(j = 0; j < K; j++) {
			erased[combinations[i][j]] = 0;
		}
		jerasure_make_decoding_matrix(K, M, W, matrix, erased, decode_matrix, dm_ids);
        print_ids(dm_ids);
		jerasure_print_matrix(decode_matrix, K, K, W);
		for(j = 0; j < K; j++) {
            all_decode[i * (K + M) + dm_ids[j]] = decode_matrix[j];
        }
		
	}
    for(i = 0; i < count; i++) {
        for(j = 0; j < K + M; j++) {
            printf("%d\t", all_decode[i * (K + M) + j]);
        }
        printf("\n");
    }
    int* all_bitmatrix = jerasure_matrix_to_bitmatrix(K + M, count ,W, all_decode);
    cout << " *************************" << endl;
    jerasure_print_bitmatrix(all_bitmatrix, count*W, (K+M)*W, W);
    cout << " *************************" << endl;

    vector<vector<int>> all_recovery;
    transform_bit_to_integer(all_bitmatrix, count*W, K + M, W, all_recovery);
    for(vector<int> tmp: all_recovery) {
        cout << tmp << endl;
    }

    vector<set<int>> xor_hash(K + M, set<int>());
    set<int> exist_hash;
    
    for(i = 0; i < K + M; i++) {
        if(all_recovery[FIRST_SELECT][i] != 0) {
            xor_hash[i].insert(all_recovery[FIRST_SELECT][i]);
        }
    }
    exist_hash.insert(FIRST_SELECT % W);
    int canDecrease = 0;
    vector<int> collect_idx;
    collect_idx.push_back(FIRST_SELECT);
    cout << "decrease: ";
    while(exist_hash.size() != W) {
        int idx;
        int decrease = getRowNum(all_recovery, xor_hash, exist_hash, W, idx);
        cout << decrease << " ";
        canDecrease += decrease;
        for(i = 0; i < K + M; i++) {
            if(all_recovery[idx][i] != 0) {
                add_xor_hash(xor_hash[i], all_recovery[idx][i]);
            }
        }
        exist_hash.insert(idx % W);
        collect_idx.push_back(idx);
    }
    cout << "**************K, M, W*************" << K << "," << M << "," << W << endl;
    cout << "origin_blocks: " << K * W << endl;
    cout << "my-tra:" << K * W - canDecrease << endl;

    cout << "collect_idx: " << collect_idx << endl;
    cout << "first_select: " << FIRST_SELECT << endl;
    for(auto xh: xor_hash) {
        for(auto x: xh) {
            cout << x << " ";
        }
        cout << endl;
    }

}


int main() {
    tan_heart_combine_reserve();
}
