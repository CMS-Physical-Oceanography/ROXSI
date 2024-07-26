%% CloseVars.m
%
% Purpose: Close all of the variable editor windows that are opened
%

desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
Titles  = desktop.getClientTitles;
for k = 1:numel(Titles)
   Client = desktop.getClient(Titles(k));     
   if ~isempty(Client) & ...
      strcmp(char(Client.getClass.getName), 'com.mathworks.mde.array.ArrayEditor')
      Client.close();
   end
end



clear Client desktop k Titles

