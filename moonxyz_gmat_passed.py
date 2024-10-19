from jplephem.spk import SPK
import numpy as np
import math

np.set_printoptions(precision=5)  #小数点后两位
kernel = SPK.open('de405.bsp')



'''
position, velocity = kernel[0,3].compute_and_differentiate(jd)  #3是地球barycenter 
print('x y z', position)    #单位km

velocity_per_second = velocity / 86400.0
print('vx vy vx',velocity_per_second)  #单位km/s
'''

'''
以下代码已与GMAT交叉验证，本代码中的儒略日对应GMAT中当日零时，坐标对应EarthMJ2000Eq下的月球坐标.
'''

def gregorian_to_julian_day(year, month, day):
    if month <= 2:
        year -= 1
        month += 12
    A = year // 100
    B = 2 - A + A // 4
    JD = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
    return JD

# 示例
year = 2023
month = 11
day = 16
jd = gregorian_to_julian_day(year, month, day)
print(f"The Julian Day for {year}-{month}-{day} is {jd}")



print()


position = kernel[3,301].compute(jd)      
position -= kernel[3,399].compute(jd) 
print('moon x y z', position)    #单位km
diyue_juli = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
print('diyue_juli', diyue_juli)
vector_A = position



print()


position = kernel[3,301].compute(jd+5)      
position -= kernel[3,399].compute(jd+5) 
print('moon x y z', position)    #单位km
diyue_juli = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
print('diyue_juli', diyue_juli)

vector_B = position




# 计算点积
dot_product = sum(a*b for a, b in zip(vector_A, vector_B))

# 计算向量的模
norm_A = math.sqrt(sum(a**2 for a in vector_A))
norm_B = math.sqrt(sum(b**2 for b in vector_B))

# 计算夹角的余弦值
cos_theta = dot_product / (norm_A * norm_B)

# 计算夹角（以弧度为单位）
theta_radians = math.acos(cos_theta)

# 将弧度转换为度
theta_degrees = math.degrees(theta_radians)

print(f"夹角的弧度: {theta_radians}")
print(f"夹角的度数: {theta_degrees}")




