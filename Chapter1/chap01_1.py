import numpy as np
from itertools import combinations

# 实体类　Solution
# 控制类　Simplex

class Solution:
    def __init__(self):
        pass

    def set_para(self, A, b, z):
        # A m*n
        # b m*1
        # z 1*(m+1)
        self.A = A
        self.b = b
        self.z = z

        self.m, self.n = A.shape
        self.x_index = [i for i in range(self.n)]

    def get_init_solution(self):
        for JB in combinations(range(self.n), self.m):
            if self._is_solution(JB):
                JB, JN = self._rearrange(JB)
                self._set_JB_JN(JB, JN)
                return True
        return False

    def _is_solution(self, JB):
        B = np.hstack([self.A[:, i] for i in JB])
        # B的行列式的值不为0，则B是一个可逆矩阵
        if np.linalg.det(B):
            return True
        return False

    def _rearrange(self, JB):
        JN = [i for i in range(self.n) if i not in JB]
        B = np.hstack([self.A[:, i] for i in JB])
        N = np.hstack([self.A[:, i] for i in JN])
        self.z = [self.z[0, i] for i in JB] + [self.z[0, i] for i in JN] + [self.z[0, -1]]
        self.z = np.matrix([self.z], dtype=float)
        # 将B转化为单位矩阵
        # 即(B|N)x=b -> (I|B'N)x=B'b
        self.A = np.hstack((np.eye(self.m), B.I * N))
        self.b = np.dot(B.I, self.b)
        # 相应的重新定义JN，JB
        self.x_index = list(JB) + list(JN)
        JB = [i for i in range(self.m)]
        JN = [i for i in range(self.m, self.n)]
        # 在目标函数中，使用非基变量替代基变量
        for i in range(self.m):
            _change = np.zeros((1, self.n + 1))
            _change[0, :self.n] = self.A[i, :]
            _change[0, -1] = -self.b[i, 0]
            self.z -= _change * self.z[0, i]
        return JB, JN

    def _set_JB_JN(self, JB, JN):
        self.JB = JB
        self.JN = JN

    def is_best(self):
        best, inf_solution = True, False
        for i in self.JN:
            sigma = self.z[0, i]
            if sigma > 0:
                best = False
            elif sigma == 0:
                inf_solution = True
        return best, inf_solution

    def get_inVar(self):
        greatest_sigma = 0
        for i in self.JN:
            sigma = self.z[0, i]
            if greatest_sigma < sigma:
                greatest_sigma = sigma
                inVar = i
        return inVar

    def get_outVar(self, inVar):
        min_ratio = self.b[self.JB[0], 0] / self.A[self.JB[0], inVar]
        outVar = self.JB[0]
        flag = False
        for i in self.JB[1:]:
            k = self.A[i, inVar]
            if k > 0:
                flag = True
                _tmp = self.b[i, 0] / k
                if _tmp < min_ratio:
                    min_ratio = _tmp
                    outVar = i
        if flag == False:
            return None
        return outVar

    def in_and_out(self, inVar, outVar):
        self.A[:, [inVar, outVar]] = self.A[:, [outVar, inVar]]
        self.x_index[outVar], self.x_index[inVar] = self.x_index[inVar], self.x_index[outVar]
        self.z[0, inVar], self.z[0, outVar] = self.z[0, outVar], self.z[0, inVar]
        B = np.hstack([self.A[:, i] for i in self.JB])
        N = np.hstack([self.A[:, i] for i in self.JN])
        # 将B转化为单位矩阵
        # 即(B|N)x=b -> (I|B'N)x=B'b
        self.A = np.hstack((np.eye(self.m), B.I * N))
        self.b = np.dot(B.I, self.b)
        # 在目标函数中，使用非基变量替代基变量

        for i in range(self.m):
            _change = np.zeros((1, self.n + 1))
            _change[0, :self.n] = self.A[i, :]
            _change[0, -1] = -self.b[i, 0]
            self.z -= _change * self.z[0, i]

    def getX(self):
        x = [0] * self.n
        for i in self.JB:
            x[self.x_index[i]] = self.b[i, 0]
        return x


class Simplex:
    def __init__(self):
        self.solution = Solution()

        # 0 正常，尚未得到最优解，继续迭代
        # 1 无解，无界解
        # 2 达到最优解
        # 3 问题有无数个最优解
        self.status = 0

    def set_para(self, A, b, z):
        # A,b,z 需以矩阵的形式输入
        self.solution.set_para(A, b, z)

    def output_result(self):
        self._main()
        if self.status == 1:
            print("此问题无界")
        elif self.status == 2:
            print("此问题有一个最优解")
        elif self.status == 3:
            print("此问题有无穷多个最优解")

    def _main(self):
        # 获得初始可行解
        self._get_init_solution()
        if self.status == 1:
            return

        while True:
            print("--------------------")
            print("z:", self.solution.z[0, -1])
            print("x:", self.solution.getX())
            # 最优性检验
            self._is_best()
            if self.status in (2, 3):
                return

            # 换入换出
            self._mainloop()
            if self.status in (1, 2):
                return

    def _get_init_solution(self):
        if self.solution.get_init_solution():
            self.status = 0
        else:
            self.status = 1

    def _is_best(self):
        best, inf_solution = self.solution.is_best()
        if best == True and inf_solution == False:
            self.status = 2
        elif best == True and inf_solution == True:
            self.status = 3
        else:
            self.status = 0

    def _mainloop(self):
        inVar = self.solution.get_inVar()
        outVar = self.solution.get_outVar(inVar)
        # 未找到换出基变量，此问题有无界解
        if outVar == None:
            self.status = 1
            return
        self.solution.in_and_out(inVar, outVar)


if __name__ == "__main__":
    s = Simplex()
    A = np.matrix([[30, 20, 1, 0, 0],
                   [5, 1, 0, 1, 0],
                   [1, 0, 0, 0, 1]])

    b = np.matrix([[160, 15, 4]]).T
    #             sigma,...,z0
    z = np.matrix([[5, 2, 0, 0, 0, 0]])
    s.set_para(A, b, z)
    s.output_result()
