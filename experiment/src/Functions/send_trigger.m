function [t_on t_off] = send_trigger(trigger, param)

t_on = GetSecs;
outp(param.io_address, trigger);

WaitSecs(param.triggerDur_inSecs);

t_off = GetSecs;
outp(param.io_address, 0)

end