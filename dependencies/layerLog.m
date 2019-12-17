function [supra,gran,infra] = layerLog(file)
%layerLog Retrieves contact information from this reportoire
%   Input BRdatafile, output three variable
switch file
    case '160418_E_mcosinteroc001'
        supra = 13:19;
        gran = 20:25;
        infra = 26:31;
    case '160420_E_mcosinteroc002'
        supra = 13:20;
        gran = 21:26;
        infra = 27:32;
    case '160421_E_mcosinteroc001'
        supra = 16:19;
        gran = 20:25;
        infra = 26:32;
    case '160422_E_mcosinteroc001'
        supra = 13:17;
        gran = 18:23;
        infra = 22:29;
    case '160423_E_mcosinteroc001'
        supra = 9:17;
        gran = 18:23;
        infra = 24:29;
    case '160423_E_mcosinteroc002'
        supra = 9:17;
        gran = 18:23;
        infra = 24:29;
    case '160425_E_mcosinteroc001'
        supra = 18:21;
        gran = 22:27;
        infra = 28:32;
    case '160429_E_mcosinteroc001'
        supra = 14:19;
        gran = 20:25;
        infra = 26:31;
    case '160429_E_mcosinteroc002'
        supra = 14:19;
        gran = 20:25;
        infra = 26:31;
    case '160502_E_mcosinteroc001'
        supra = 13:20;
        gran = 21:26;
        infra = 27:32;
    case '160505_E_mcosinteroc001'
        supra = 14:22;
        gran = 23:28;
        infra = 29:32;
    case '160505_E_mcosinteroc002'
        supra = 14:22;
        gran = 23:28;
        infra = 29:32;
    case '160510_E_mcosinteroc001'
        supra = 19:23;
        gran = 24:29;
        infra = 30:32;
    case '160510_E_mcosinteroc002'
        supra = 19:23;
        gran = 24:29;
        infra = 30:32;
    case '160512_E_mcosinteroc001'
        supra = 10:16;
        gran = 17:22;
        infra = 23:32;
    case '160512_E_mcosinteroc002'
        supra = 10:16;
        gran = 17:22;
        infra = 23:32;
    case '160523_E_mcosinteroc002'
        supra = 10:13;
        gran = 14:19;
        infra = 20:24;
    case '170719_I_mcosinteroc003'
        supra = 12:14;
        gran = 15:20;
        infra = 21:32;
    case '170724_I_mcosinteroc001'
        supra = 5:14;
        gran = 15:20;
        infra = 21:32;
end

