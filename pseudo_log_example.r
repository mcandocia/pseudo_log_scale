source('pseudo_log_scale.r')

set.seed(2020 * 2019)
n_subjects = 8
n_time = 100

df_base = expand.grid(
  subject_id = 1:n_subjects,
  time = 1:n_time
)

time_comp = data.frame(
  time = 1:n_time,
  time_value = 1+cumsum(rnorm(n_time)/30)
)

subject_comp = data.frame(
  subject_id = 1:n_subjects,
  subject_value = 10^seq(2,6, length.out=n_subjects)
)

df = df_base %>%
  inner_join(
    time_comp
  ) %>%
  inner_join(
    subject_comp
  ) %>%
  mutate(
    value = (time_value + rnorm(n())/100) * subject_value + rnorm(n()),
    subject_id = as.character(subject_id)
  ) %>%
  group_by(
    subject_id
  ) %>%
  mutate(
    value_diff = c(NA_real_, diff(value))
  ) %>%
  ungroup()

ggplot(df %>% filter(time <= 20)) + 
  geom_line(
    aes(x=time, y=value_diff, color=subject_id)
  )

ggplot(df %>% filter(time <= 20)) + 
  geom_line(
    aes(x=time, y=value_diff, color=subject_id)
  ) + 
  scale_y_continuous(
    trans=make_lal_trans('lal',1,10)
  )

