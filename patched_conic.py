from jplephem.spk import SPK
import numpy as np
import math
from constant import *



np.set_printoptions(precision=5)  #小数点后两位
kernel = SPK.open('de405.bsp')



'''
position, velocity = kernel[0,3].compute_and_differentiate(jd)  #3是地球barycenter 
print('x y z', position)    #单位km

velocity_per_second = velocity / 86400.0
print('vx vy vx',velocity_per_second)  #单位km/s
'''



def gregorian_to_julian_day(year, month, day):
    if month <= 2:
        year -= 1
        month += 12
    A = year // 100
    B = 2 - A + A // 4
    JD = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day + B - 1524.5
    return JD



print()

def patched_conic(jd, flytime, y_value, ramda2, r_s, leo):
  position = kernel[3,301].compute(jd)      
  position -= kernel[3,399].compute(jd) 
  #print('moon x y z', position)    #单位km
  diyue_juli = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
  #print('diyue_juli', diyue_juli)
  vector_A = position
  #print()
  position = kernel[3,301].compute(jd+5)   #5d后moon位置    
  position -= kernel[3,399].compute(jd+5) 
  #print('moon x y z', position)    #单位km
  diyue_juli = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
  #print('diyue_juli', diyue_juli)
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
  #print(f"夹角的弧度: {theta_radians}")
  #print(f"夹角的度数: {theta_degrees}")
  # 计算单位向量
  unit_vector_b = vector_B / np.linalg.norm(vector_B)
  r_park = r_earth + leo             
  r_e0 = - r_park * unit_vector_b                  
  #print('r_e,0',r_e0)
  '''
  r_e0_gaodu = np.sqrt(r_e0[0]**2 + r_e0[1]**2 + r_e0[2]**2)
  print('r_e0_gaodu', r_e0_gaodu)

  '''

  # 将 vb 归一化作为 x 轴
  va = vector_A / np.linalg.norm(vector_A)
  vb = vector_B / np.linalg.norm(vector_B)
  x_axis = vb
  # 计算法向量（z轴）
  z_axis = np.cross(va, vb)
  z_axis = z_axis / np.linalg.norm(z_axis)  # 确保z轴也是单位向量
  # 计算y轴
  y_axis = np.cross(z_axis, x_axis)    #x轴向右，y轴向上
  y_axis = y_axis / np.linalg.norm(y_axis)  # 确保y轴也是单位向量
  # 确保坐标轴正交
  assert np.allclose(np.dot(x_axis, y_axis), 0)
  assert np.allclose(np.dot(y_axis, z_axis), 0)
  assert np.allclose(np.dot(z_axis, x_axis), 0)

  # 构造向量 [0, y, 0]
  
  vector_in_plane = y_value * y_axis   
  #print('v_e,0',vector_in_plane)
  '''
  x_value = - r_park
  vector2_in_plane = x_value * x_axis
  print(vector2_in_plane)
  '''
  vector_B_magnitude = np.linalg.norm(vector_B)
  vector_B_in_plane_coordinates = np.array([vector_B_magnitude, 0, 0])
  #print(vector_B_in_plane_coordinates)
  r_e0_magnitude = np.linalg.norm(r_e0)
  r_e0_in_plane_coordinates = np.array([-r_e0_magnitude, 0, 0])    #书里的r_e,0
  #print(r_e0_in_plane_coordinates)
  v_e1 = np.array([0, y_value, 0])       #书里的v_e,1
  #print(v_e1)
  v_e1_magnitude = np.linalg.norm(v_e1)
  #print(v_e1_magnitude)

  eng_e = v_e1_magnitude**2/2 - mu_earth/r_e0_magnitude
  #print(eng_e)
  a_e = - mu_earth/(2*eng_e)
  #print(a_e)
  phi_e1 = 0
  h_e = r_e0_magnitude*v_e1_magnitude*math.cos(phi_e1)
  #print(h_e)
  p_e = h_e**2/mu_earth
  #print('p_e',p_e)
  e_e = math.sqrt(1-p_e/a_e)
  #print('e_e',e_e)
  
  position = kernel[3,301].compute(jd+flytime)    
  position -= kernel[3,399].compute(jd+flytime) 
  #print('moon x y z', position)    #单位km
  diyue_juli = np.sqrt(position[0]**2 + position[1]**2 + position[2]**2)
  #print('diyue_juli', diyue_juli)
  vector_C = position
  # 计算vc在x_axis和y_axis上的投影
  vc_x = np.dot(vector_C, x_axis)
  vc_y = np.dot(vector_C, y_axis)
  vc_z = np.dot(vector_C, z_axis)
  # vc 在新坐标系下的坐标
  vc_coords = np.array([vc_x, vc_y, vc_z])
  vc_magnitude = np.linalg.norm(vc_coords)
  #print(vc_magnitude)
  
  ramda_2 = ramda2 * math.pi / 180
  r_e2 = math.sqrt(vc_magnitude**2 + r_s**2 - 2*vc_magnitude*r_s*math.cos(ramda_2))
  #print(r_e2)


  v_e2 = math.sqrt(2*(eng_e+mu_earth/r_e2))
  #print('v_e2',v_e2)
  phi_e2 = math.acos(h_e/(r_e2*v_e2))
  #print(phi_e2 * (180 / math.pi))
  qq = (1 - e_e) / (1 + e_e)          #中间量，不用揣测代号含义
  # 计算f1
  f1 = math.acos((p_e - r_e0_magnitude) / (e_e * r_e0_magnitude))
  #print()
  # 计算E1
  E1 = 2 * math.atan(math.sqrt(qq) * math.tan(f1 / 2))
  # 计算f2
  f2 = math.acos((p_e - r_e2) / (e_e * r_e2))
  # 计算E2
  E2 = 2 * math.atan(math.sqrt(qq) * math.tan(f2 / 2))
  #print(f1)
  #print(f2)
  #print(E1)
  #print(E2)
  # 计算t_e
  t_e = ((E2 - E1) - e_e * math.sin(E2 - E1)) * math.sqrt(a_e**3/mu_earth)    #妈的书上写错了，艹
  #print(t_e / 3600)         #单位是s
  v_min = math.sqrt(2*mu_earth*r_e2/(r_e0_magnitude*(r_e0_magnitude+r_e2)))
  #print(v_min)
  psi2 = math.asin(math.sin(ramda_2)*r_s/r_e2)
  position1, velocity1 = kernel[3,301].compute_and_differentiate(jd+flytime)  #3是地球barycenter 
  position2, velocity2 = kernel[3,399].compute_and_differentiate(jd+flytime)
  v_moon = np.linalg.norm((velocity1-velocity2) / 86400.0)
  #print('v_moon',v_moon)  #单位km/s
  v_m2 = math.sqrt(v_e2**2 + v_moon**2 -2*v_e2*v_moon*math.cos(phi_e2-psi2))
  #print(v_m2)
  qqq = (v_moon*math.cos(ramda_2)-v_e2*math.cos(ramda_2+psi2-phi_e2))/v_m2    #中间量，不用揣测代号含义
  epsilon2 = math.asin(qqq)
  #print(epsilon2* (180 / math.pi))


  eng_moon = v_m2**2/2-mu_moon/r_s
  a_moon = -mu_moon/(2*eng_moon)
  h_moon = r_s*v_m2*math.sin(epsilon2)
  p_moon = h_moon**2/mu_moon
  e_moon = math.sqrt(1-p_moon/a_moon)
  r_moonpe = p_moon/(1+e_moon)
  v_moonpe = math.sqrt(2*(eng_moon+mu_moon/r_moonpe))

  '''
  print()
  print(eng_moon)
  print(a_moon)
  
  print(p_moon)
  
  '''
  print('r_moonpe',r_moonpe, 'km')
  print('v_moonpe',v_moonpe, 'km/s')
  print('e_moon',e_moon)
  print('h_moon',h_moon)
  return r_moonpe,e_moon,h_moon,t_e
  '''
  在月球影响球内就是双曲线轨道了，需要用双曲线轨道的偏近点角的相关公式计算飞行时间
  '''

def find_jiao(r_moonpe, r_moon, mu_moon, e_moon, r_s):
    jiao_found = None

    # 遍历角度jiao从1到100
    for jiao in range(-300, 0):
        # 计算r_loop
        r_loop = (r_moonpe + r_moon)**2 / (mu_moon * (1 + e_moon * math.cos(math.radians(jiao))))
        
        # 检查条件是否满足
        if abs(r_loop) - r_s < 0.1:
            jiao_found = jiao
            break  # 找到满足条件的jiao后退出循环

    # 循环结束后检查是否找到了满足条件的jiao
    if jiao_found is not None:
        return jiao_found
    else:
        # 如果没有找到满足条件的jiao，可以返回一个提示或者None
        return None


# 示例
year = 2020
month = 12
day = 1
jd = gregorian_to_julian_day(year, month, day)
print(f"The Julian Day for {year}-{month}-{day} is {jd}")


flytime = 4.67    #抵达从jd起算flyime天后的月球
y_value = -10.933  # 加速后的速度

print('发射速度（符号仅代表方向）',y_value,'km/s')
r_s = 66200       #km,月球影响球半径
ramda2 = -30    #degree
leo = 200        #表示200km的LEO轨道上
r_moonpe,e_moon,h_moon,t_e = patched_conic(jd, flytime, y_value, ramda2, r_s, leo)

# 你可以这样调用这个函数：
result = find_jiao(r_moonpe, r_moon, mu_moon, e_moon, r_s)
#print("Found jiao:", result)
qiu = math.sqrt((e_moon-1)/(e_moon+1))*math.tan(math.radians(result))   #中间量，不用揣测代号含义
#print(qiu)
F_half = math.atanh(qiu)
F = 2 * F_half
M_h = e_moon*F-F
shijian = h_moon**3 * M_h/(mu_moon**2 * (e_moon**2 - 1)**1.5)
#print(-shijian/3600)

total_flytime = -shijian/3600 +t_e / 3600
print('total flytime',total_flytime,'hours')  #hour
fashe_shike = flytime - total_flytime/24
print('发射时刻为',jd+fashe_shike)


