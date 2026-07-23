# External Int2: прямой аналитический разбор

## 1. Исходный интеграл

Пусть

\[
r=\frac{s}{t}>0,
\]

и чистая параметрическая часть равна

\[
\mathcal J_2(\varepsilon,r)=
\int_{\mathbb R_+^3}
 x_2^{1+\varepsilon}(1+x_2)^\varepsilon(1+x_5)^\varepsilon
 (1+x_7)^{-1-\varepsilon}
 (1+x_7+x_2x_7+r x_2x_5)^{-1+\varepsilon}
 \,dx_2dx_5dx_7.
\]

Полный объект:

\[
I_2=P_2(\varepsilon,t)\,\mathcal J_2(\varepsilon,r),
\]

\[
P_2=e^{2\varepsilon\gamma_E}t^{-3-\varepsilon}
\frac{\Gamma(1-\varepsilon)\Gamma(-\varepsilon)^3\Gamma(\varepsilon)}
{\Gamma(-1-3\varepsilon)\Gamma(-2\varepsilon)}.
\]

## 2. Точное интегрирование по x7

Для

\[
A=1+r x_2x_5,\qquad B=1+x_2
\]

имеем

\[
\int_0^\infty(1+x_7)^{-1-\varepsilon}(A+Bx_7)^{-1+\varepsilon}dx_7
=\frac{B^\varepsilon-A^\varepsilon}{\varepsilon(B-A)}.
\]

Поэтому

\[
\mathcal J_2=
\frac1\varepsilon\int_0^\infty\!dx_2dx_5\,
 x_2^\varepsilon(1+x_2)^\varepsilon(1+x_5)^\varepsilon
 \frac{(1+x_2)^\varepsilon-(1+r x_2x_5)^\varepsilon}{1-rx_5}.
\]

После интегрирования по x2:

\[
\mathcal J_2=\frac{G_1(\varepsilon)}{\varepsilon}Q(\varepsilon,r),
\]

\[
G_1(\varepsilon)=
\frac{\Gamma(1+\varepsilon)\Gamma(-1-3\varepsilon)}{\Gamma(-2\varepsilon)},
\]

\[
Q(\varepsilon,r)=\int_0^\infty
(1+x)^\varepsilon
\frac{1-{}_2F_1(-\varepsilon,1+\varepsilon;-2\varepsilon;1-rx)}{1-rx}
\,dx.
\]

## 3. Точное дифференциальное уравнение

Обозначим

\[
C(\varepsilon)=
\frac{\Gamma(-2\varepsilon)\Gamma(-1-2\varepsilon)}
{\Gamma(-\varepsilon)\Gamma(-1-3\varepsilon)}.
\]

Из гипергеометрического уравнения, с учетом обоих граничных вкладов
z=1 и z=-\infty, получается

\[
\begin{aligned}
r^2(1+r)^2 Q''
&+r(1+r)\bigl(3(1+r)+\varepsilon(3-r)\bigr)Q'\\
&+\bigl((1+r)^2+\varepsilon(3+r-r^2)+\varepsilon^2(2-r)\bigr)Q\\
&=\varepsilon r-(1+2\varepsilon)
\left(1-C(\varepsilon)+C(\varepsilon)r^{-\varepsilon}\right).
\end{aligned}
\]

Это уравнение численно проверено непосредственно на одномерном интеграле Q.

Положим S=rQ и

\[
S=\frac1\varepsilon+s_0+\varepsilon s_1+\varepsilon^2s_2+\varepsilon^3s_3+\cdots.
\]

Рекурсия имеет вид

\[
\mathcal D_0s_n=R_n-\mathcal D_1s_{n-1}-\mathcal D_2s_{n-2},
\]

\[
\mathcal D_0=(1+r)^2(r\partial_r)',\qquad
\mathcal D_1=(1+r)(3-r)\partial_r-1,\qquad
\mathcal D_2=\frac2r-1.
\]

## 4. Решение рекурсии

Пусть L=log r и

\[
H_{a_1,\ldots,a_n}(r)=G(a_1,\ldots,a_n;r),\qquad a_i\in\{0,-1\}.
\]

Тогда

\[
s_0=-\frac34L,
\]

\[
s_1=\frac18L^2-\frac{\pi^2}{4},
\]

\[
s_2=\frac18L^3-\frac12H_{-1,0,0}
+\frac{11\pi^2}{24}L-\frac{\pi^2}{4}H_{-1}-\frac12\zeta_3,
\]

\[
\begin{aligned}
s_3={}&-\frac{11}{96}L^4-\frac{7\pi^2}{16}L^2+L\zeta_3
+\pi^2H_{0,-1}+2H_{0,-1,0,0}\\
&+\frac{\pi^2}{3}H_{-1,0}+\frac32H_{-1,0,0,0}
-\frac{\pi^2}{4}H_{-1,-1}-\frac12H_{-1,-1,0,0}\\
&-\frac12H_{-1}\zeta_3-\frac{\pi^4}{9}.
\end{aligned}
\]

Подстановка этих функций в рекурсию дает точный ноль по каждому порядку.
Высокоточная одномерная квадратура независимо воспроизводит константы интегрирования.

## 5. Полный Laurent-ряд через eps^0

Пусть

\[
T=\log t,\qquad L=\log\frac{s}{t},\qquad r=\frac{s}{t}.
\]

Тогда

\[
I_2=\frac1{s t^2}
\left[
-\frac4{\varepsilon^4}
+\frac{F_{-3}}{\varepsilon^3}
+\frac{F_{-2}}{\varepsilon^2}
+\frac{F_{-1}}{\varepsilon}
+F_0+O(\varepsilon)
\right],
\]

где

\[
F_{-3}=3L+4T,
\]

\[
F_{-2}=\frac{5\pi^2}{3}-\frac12L^2-3LT-2T^2,
\]

\[
\begin{aligned}
F_{-1}={}&2H_{-1,0,0}-\frac12L^3+\frac12L^2T
+\frac32LT^2+\frac23T^3\\
&-\frac{7\pi^2}{3}L-\frac{5\pi^2}{3}T
+\pi^2H_{-1}+\frac{62}{3}\zeta_3,
\end{aligned}
\]

и

\[
\begin{aligned}
F_0={}&\frac{11}{24}L^4+\frac12L^3T
+L^2\left(-\frac14T^2+\frac{11\pi^2}{6}\right)\\
&+L\left(-\frac12T^3+\frac{7\pi^2}{3}T-18\zeta_3\right)
-\frac16T^4+\frac{5\pi^2}{6}T^2\\
&-2T H_{-1,0,0}-\frac{62}{3}T\zeta_3
-4\pi^2H_{0,-1}-8H_{0,-1,0,0}\\
&-\frac{4\pi^2}{3}H_{-1,0}-6H_{-1,0,0,0}
+\pi^2H_{-1,-1}+2H_{-1,-1,0,0}\\
&+H_{-1}(-\pi^2T+2\zeta_3)+\frac{23\pi^4}{45}.
\end{aligned}
\]

Здесь все H-функции имеют аргумент r=s/t.

После обратного масштабирования GPL

\[
G(-t,0,0;s)=H_{-1,0,0}(s/t),\qquad
G(0,-t,0,0;s)=H_{0,-1,0,0}(s/t),
\]

и использования

\[
G(0,\ldots,0;z)=\frac{\log^n z}{n!},
\]

полученный ряд совпадает с AnsvInt2 из notebook по всем порядкам
\varepsilon^{-4},\ldots,\varepsilon^0.

## 6. Независимые проверки

- точное интегрирование по x7;
- численная проверка гипергеометрических connection formulas;
- точное дифференциальное уравнение для Q;
- точная проверка рекурсии для s0,...,s3;
- прямая одномерная квадратура Q при нескольких eps и r;
- численное сравнение компактной формулы с исходным GPL-выражением notebook при произвольных положительных s,t.

## Вывод

Int2 аналитически вычисляется без обращения к найденному non-LF IBP-базису:
точное интегрирование по x7 сводит задачу к одномерному гипергеометрическому интегралу,
для которого существует замкнутое дифференциальное уравнение по r=s/t.
Его epsilon-рекурсия порождает GPL/HPL-алфавит {0,-1} и воспроизводит весь предоставленный Laurent-ряд через конечную часть.
