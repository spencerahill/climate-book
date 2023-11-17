#!/usr/bin/env python
# coding: utf-8

# # Introduction to feedbacks

# ____________
# <a id='section3'></a>
# 
# ## 3. The feedback concept
# ____________
# 
# A concept borrowed from electrical engineering. You have all heard or used the term before, but we’ll try take a more precise approach today.
# 
# A feedback occurs when a portion of the output from the action of a system is added to the input and subsequently alters the output:
# 
# 
# <img src="http://www.atmos.albany.edu/facstaff/brose/classes/ENV415_Spring2018/images/feedback_sketch.png" width="500">
# 
# The result of a loop system can either be **amplification** or **dampening** of the process, depending on the sign of the gain in the loop, which we will denote $f$.
# 
# We will call amplifying feedbacks **positive** ($f>0$) and damping feedbacks **negative** ($f<0$).
# 
# We can think of the “process” here as the entire climate system, which contains many examples of both positive and negative feedback.

# ### Example: the water vapor feedback
# 
# The capacity of the atmosphere to hold water vapor (saturation specific humidity) increases exponentially with temperature. Warming is thus accompanied by moistening (more water vapor), which leads to more warming due to the enhanced water vapor greenhouse effect.
# 
# **Positive or negative feedback?**
# 
# 
# ### Example: the ice-albedo feedback
# 
# Colder temperatures lead to expansion of the areas covered by ice and snow, which tend to be more reflective than water and vegetation. This causes a reduction in the absorbed solar radiation, which leads to more cooling. 
# 
# **Positive or negative feedback?**
# 
# *Make sure it’s clear that the sign of the feedback is the same whether we are talking about warming or cooling.*
# 
# 

# _____________
# <a id='section4'></a>
# ## 4. Climate feedback: some definitions
# ____________
# 
# We start with an initial radiative forcing , and get a response
# $$ \Delta T_0 = \frac{\Delta R}{\lambda_0} $$
# 
# 
# Now consider what happens in the presence of a feedback process. For a concrete example, let’s take the **water vapor feedback**. For every degree of warming, there is an additional increase in the greenhouse effect, and thus additional energy added to the system.
# 
# Let’s denote this extra energy as 
# $$ f \lambda_0 \Delta T_0 $$
# 
# where $f$ is the **feedback amount**, a number that represents what fraction of the output gets added back to the input. $f$ must be between $-\infty$ and +1. 
# 
# For the example of the water vapor feedback, $f$ is positive (between 0 and +1) – the process adds extra energy to the original radiative forcing.
# 
# The amount of energy in the full "input" is now
# 
# $$ \Delta R + f \lambda_0 \Delta T_0 $$
# 
# or
# 
# $$ (1+f) \lambda_0 \Delta T_0 $$
# 
# But now we need to consider the next loop. A fraction $f$ of the additional energy is also added to the input, giving us
# 
# $$ (1+f+f^2) \lambda_0 \Delta T_0 $$
# 
# 

# and we can go round and round, leading to the infinite series
# 
# $$ (1+f+f^2+f^3+ ...) \lambda_0 \Delta T_0 = \lambda_0 \Delta T_0 \sum_{n=0}^{\infty} f^n $$
# 
# Question: what happens if $f=1$?
# 
# It so happens that this infinite series has an exact solution
# 
# $$ \sum_{n=0}^{\infty} f^n = \frac{1}{1-f} $$
# 
# So the full response including all the effects of the feedback is actually
# 
# $$ \Delta T = \frac{1}{1-f} \Delta T_0 $$
# 
# This is also sometimes written as 
# $$ \Delta T = g \Delta T_0 $$
# 
# where 
# 
# $$ g = \frac{1}{1-f} = \frac{\Delta T}{\Delta T_0} $$
# 
# is called the **system gain** -- the ratio of the actual warming (including all feedbacks) to the warming we would have in the absence of feedbacks.
# 
# So if the overall feedback is positive, then $f>0$ and $g>1$.
# 
# And if the overall feedback is negative?
# 
# 

# In[1]:


# code to plot feedback...


# In[ ]:




